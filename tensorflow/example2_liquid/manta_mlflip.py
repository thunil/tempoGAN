# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Fluid implicit particle (FLIP) with Machine Learning
#
# ----------------------------------------------------------------------------

import os, sys, argparse, pickle
assertNumpy() # make sure mantaflow is compiled with the NUMPY option, ie, "cmake ... -DNUNPY=1"

def path_to_frame(outdir, frame):
	return '{}/{:05d}/'.format(outdir, frame)

def save_frame(outdir, frame, saving_funcs):
	if (outdir is None) or (params['frame_saved']==frame): return

	path = path_to_frame(outdir, frame)
	os.path.isdir(path) or os.makedirs(path)
	for save, name in saving_funcs: save(path+name, notiming=True)

	params['frame_saved'] = frame
	print('Frame #{} was saved in {}\n'.format(frame, path))

parser = argparse.ArgumentParser(description='FLIP with ML', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(      '--load',   default='../data/mlflip-tf/',  help='path to the trained tensorflow model directory')
parser.add_argument('-w', '--window', default=1, type=int, help='window size for sampling features; 1 (default) means 3x3, 2 means 5x5, so on.')
parser.add_argument(      '--t_end', default=6.0, type=float, help='end time of simulation.')
pargs = parser.parse_known_args()[0]

onphi  = True
ongeom = False

if pargs.load is None: sys.exit('You have to specify the path to the trained model.')
pargs.load = os.path.normpath(pargs.load)

with open(pargs.load + '/run_args.pickle', 'rb') as f: tfopt = pickle.load(f)
with open(pargs.load + '/scale.pickle', 'rb') as f: scale = pickle.load(f)

import numpy as np
dtype_real = np.float32         # NOTE: if double precision, use float64
dtype_int  = np.int32           # NOTE: if int in C is 64bits, use int64
np.random.seed(1)

import tensorflow as tf
tf_sess = tf.InteractiveSession()

import tf_network
dlayers = list(map(int, tfopt['dnet'].split('-')))
mlayers = list(map(int, tfopt['mnet'].split('-')))
dact    = list(map(tf_network.parse_act, tfopt['dact'].split('-')))
mact    = list(map(tf_network.parse_act, tfopt['mact'].split('-')))
x       = tf.placeholder(tf.float32, shape=[None, dlayers[0]], name='x-input')
y_,  y  = tf_network.build_network(dlayers, dact, input_x_holder=x, bn=tfopt['bn'], is_training=False, scope='detector/')[1:]
y2_, y2 = tf_network.build_network(mlayers, mact, input_x_holder=x, bn=tfopt['bn'], is_training=False, scope='modifier/')[1:]
if tfopt['mve']:
	sd  = tf_network.build_network(mlayers, mact, input_x_holder=x, input_y_holder=y2_, bn=tfopt['bn'], is_training=False, scope='modifier_var/')[2]

tf_saver = tf.train.Saver()
modelfile = pargs.load + '/model.ckpt'
tf_saver.restore(tf_sess, modelfile)
print('Pre-trained model {} loaded\n'.format(modelfile))

import manta as mt
mt.tFluid    = FlagFluid
mt.tObstacle = FlagObstacle

nogui       = False
pause       = False
output      = None
#output      = '../data/manta-mlflip'
savingFuncs = []

# default solver parameters
params                = {}
params['dim']         = 2                  # dimension
params['sres']        = 2                  # particle sampling resolution per cell
params['dx']          = 1.0/params['sres'] # particle spacing (= 2 x radius)
params['res']         = 75                 # reference resolution
params['len']         = 1.0                # reference length
params['bnd']         = 2                  # boundary cells
params['grav']        = 0                  # applied gravity (mantaflow scale); recomputed later
params['gref']        = -9.8               # real-world gravity
params['jitter']      = 0.2                # jittering particles
params['cgaccuracy']  = 1e-3               # cg solver's threshold
params['fps']         = 24
params['t_end']       = pargs.t_end        # default 6.0
params['sdt']         = None
params['frame_saved'] = -1
params['frame_last']  = -1

scaleToManta = float(params['res'])/params['len']
params['gs']    = [ int(params['res']*1.5+params['bnd']*2), int(params['res']+params['bnd']*2), int(params['res']+params['bnd']*2 if params['dim']==3 else 1) ]
params['grav']  = params['gref']*scaleToManta

#params['stref'] = 0.073 # surface tension (reference scale [m]; e.g., 0.073)
#params['stens'] = params['stref']*scaleToManta

s             = Solver(name='MLFLIP', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

# prepare grids and particles
gFlags   = s.create(FlagGrid)
_dummyV_ = s.create(MACGrid) # hide velocities for UI...
_dummyR_ = s.create(RealGrid) # hide pressure for UI...

gV       = s.create(MACGrid)
gVold    = s.create(MACGrid)
gP       = s.create(RealGrid)
gPhi     = s.create(LevelsetGrid)
gFlagTmp = s.create(FlagGrid)

pp      = s.create(BasicParticleSystem)
gIdxSys = s.create(ParticleIndexSystem)
gIdx    = s.create(IntGrid)

pT    = pp.create(PdataInt)    # main particle flags
pV    = pp.create(PdataVec3)   # particle velocity
pVtmp = pp.create(PdataVec3)
pVtm2 = pp.create(PdataVec3)
pVtm3 = pp.create(PdataVec3)
pItmp = pp.create(PdataInt)   # typically used for temporary flags

mesh = s.create(name='mesh', type=Mesh) if (params['dim']==3 and not nogui) else None

savingFuncs.append([pp.save, 'particles.uni'])
savingFuncs.append([pV.save, 'particlesVel.uni'])
savingFuncs.append([pT.save, 'particlesType.uni'])

fv_N_axis = 2*pargs.window + 1
fv_N_stn  = fv_N_axis*fv_N_axis*(fv_N_axis if params['dim']==3 else 1)
fv_N_row  = params['dim']*fv_N_stn + (fv_N_stn if onphi else 0) + (fv_N_stn if ongeom else 0)
fv_vscale = params['len']/float(params['res'])

# boundary
gFlags.initDomain(params['bnd']-1)

# fluid dam setup
a = vec3(params['res']*0.2+params['bnd'], params['res']*0.2+params['bnd'], params['res']*0.2+params['bnd'] if (params['dim']==3) else 0)
b = vec3(params['res']*0.2, params['res']*0.2, params['res']*0.2 if (params['dim']==3) else params['gs'][2])
fld = s.create(Box, center=a, size=b)

begin = pp.pySize()
sampleShapeWithParticles(shape=fld, flags=gFlags, parts=pp, discretization=params['sres'], randomness=params['jitter'], notiming=True)
end = pp.pySize()
pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)
markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle)

gui = None
if not nogui:
	gui = Gui()
	gui.show()
	gui.setCamPos(0,0.1,-1.2)
	#gui.toggleHideGrids()
	if pause: gui.pause()

if output:
	save_frame(output, s.frame, savingFuncs)
	with open(output+'/params.pickle', 'wb') as f: pickle.dump(params, f)    

stats = {'candidate': 0, 'decision': 0, 'reverted': 0, 'splashed': 0}
np_pTimer = np.zeros(pp.pySize(), dtype=dtype_real)
while (s.timeTotal<params['t_end']): # main loop

	mapPartsToMAC(vel=gV, flags=gFlags, velOld=gVold, parts=pp, partVel=pV, ptype=pT, exclude=FlagEmpty)
	if params['sdt'] is None: s.adaptTimestep(gV.getMax())
	else: s.adaptTimestepByDt(params['sdt'])

	addGravityNoScale(flags=gFlags, vel=gV, gravity=vec3(0, params['grav'], 0))

	gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
	unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)

	setWallBcs(flags=gFlags, vel=gV)
	solvePressure(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'], phi=gPhi)
	setWallBcs(flags=gFlags, vel=gV)
	extrapolateMACSimple(flags=gFlags, vel=gV)

	# BEGIN: machine learning part >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	# get candidate particles
	gFlagTmp.copyFrom(gFlags)
	extendRegion(flags=gFlagTmp, region=FlagEmpty, exclude=FlagObstacle, depth=1)
	pItmp.copyFrom(pT)
	setPartType(parts=pp, ptype=pItmp, mark=0, stype=FlagEmpty, flags=gFlagTmp, cflag=FlagEmpty|FlagFluid) # already individual? then, kill
	setPartType(parts=pp, ptype=pItmp, mark=FlagEmpty, stype=FlagFluid, flags=gFlagTmp, cflag=FlagEmpty)   # mark surface particles
	candidate = np.zeros(pp.pySize(), dtype=dtype_int)
	copyPdataToArrayInt(target=candidate, source=pItmp)

	candidate = (candidate==FlagEmpty) # turn into bool array
	N_candidate = np.count_nonzero(candidate)
	stats['candidate'] += N_candidate

	# extract features -> numpy array
	inputs_c = np.zeros(pp.pySize()*fv_N_row, dtype=dtype_real)
	off_feature = 0
	extractFeatureVel(fv=inputs_c, N_row=fv_N_row, off_begin=off_feature, p=pp, vel=gV, scale=fv_vscale, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
	off_feature = params['dim']*fv_N_stn
	if onphi:
		extractFeaturePhi(fv=inputs_c, N_row=fv_N_row, off_begin=off_feature, p=pp, phi=gPhi, scale=1.0, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
		off_feature += fv_N_stn
	if ongeom:
		extractFeatureGeo(fv=inputs_c, N_row=fv_N_row, off_begin=params['dim']*fv_N_stn, p=pp, flag=gFlags, scale=1.0, ptype=pItmp, exclude=FlagObstacle|FlagFluid, window=pargs.window)
		off_feature += fv_N_stn

	inputs_c = inputs_c.reshape((-1, fv_N_row))[candidate]

	# approximate using tf: detection and modification
	if tfopt['mve']:    dtct_c, dv_c, appx_s_c = tf_sess.run([y, y2, sd], feed_dict={x: inputs_c})
	else:               dtct_c, dv_c           = tf_sess.run([y, y2],     feed_dict={x: inputs_c})

	# classify splashes, compute bool array from probabilities
	# when 0th value is larger, it means splashing
	dtct_bool = (np.argmax(dtct_c, axis=1)==0) 
	#dtct_bool = np.zeros(dtct_bool.size, dtype=dtype_int) # uncomment to turn the splash model off (for testing)

	# mark new decision for the splashing particles (as FlagEmpty)
	N_splashing = np.count_nonzero(dtct_bool)
	decision = np.zeros(pp.pySize(), dtype=dtype_int)
	# dtct_bool only contains info for surface particles, map back to full particle system, and convert bool to manta flag
	decision[candidate] = dtct_bool
	decision[(decision==1)] = FlagEmpty
	decision[(decision!=FlagEmpty)] = FlagFluid
	stats['decision'] += N_splashing

	# calculate velocity modification
	if params['dim']==2:
		dv_c = np.append(dv_c, np.zeros((N_candidate, 1), dtype=dtype_real), axis=1)
		if tfopt['mve']: appx_s_c = np.append(appx_s_c, np.zeros((N_candidate, 1), dtype=dtype_real), axis=1)

	if tfopt['mve']: dv_c += appx_s_c*np.random.normal(size=(N_candidate, 3))

	dv = np.zeros((pp.pySize(), 3), dtype=dtype_real)
	dv[candidate] = dv_c * scale['modvel']

	# END: machine learning part <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	# update velocity; general update from FLIP
	flipVelocityUpdate(vel=gV, velOld=gVold, flags=gFlags, parts=pp, partVel=pV, flipRatio=0.97, ptype=pT, exclude=FlagObstacle|FlagEmpty)

	# update velocity; individual update for Lagrangian particles
	addForcePvel(vel=pV, a=vec3(0, params['grav'], 0), dt=s.timestep, ptype=pT, exclude=FlagObstacle|FlagFluid)

	# 0. let's predict splash in framestep, which is the same with the training data's timestep
	pp.getPosPdata(target=pVtmp)
	c_dt = s.timestep
	s.timestep = s.frameLength

	# 1. mark fluid region after moving without splash-correction
	gFlagTmp.copyFrom(gFlags)
	pp.advectInGrid(flags=gFlagTmp, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagObstacle|FlagFluid)
	markFluidCells(parts=pp, flags=gFlagTmp, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	extendRegion(flags=gFlagTmp, region=FlagEmpty, exclude=FlagObstacle, depth=1)

	# 2. try to move splashing particles only, so check if it's really splashing; revert the wrong decisions
	pp.setPosPdata(source=pVtmp)
	pVtm2.copyFrom(pV)
	copyArrayToPdataVec3(target=pVtm3, source=dv.reshape(-1, 1))
	pVtm2.add(pVtm3)
	copyArrayToPdataInt(target=pItmp, source=decision)
	eulerStep(parts=pp, vel=pVtm2, ptype=pItmp, exclude=FlagObstacle|FlagFluid)
	setPartType(parts=pp, ptype=pItmp, mark=FlagFluid, stype=FlagEmpty, flags=gFlagTmp, cflag=FlagFluid|FlagObstacle) # empty -> fluid if they are not acturally splashing.
	copyPdataToArrayInt(target=decision, source=pItmp)

	# 3. final decision and velocity modification
	stats['splashed'] += np.count_nonzero(decision==FlagEmpty)
	dv[(decision!=FlagEmpty)] = 0
	np_pTimer[(decision==FlagEmpty)] = s.frameLength # set judgement timer

	# 4. roll-back
	s.timestep = c_dt
	pp.setPosPdata(source=pVtmp)

	# mark splashing particles and modify the velocities so that they can flow individually
	copyArrayToPdataInt(target=pItmp, source=decision)
	pT.setConstIntFlag(s=FlagEmpty, t=pItmp, flag=FlagEmpty)
	copyArrayToPdataVec3(target=pVtm2, source=dv.reshape(-1, 1))
	pV.add(pVtm2)

	# update position
	pp.advectInGrid(flags=gFlags, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagObstacle|FlagFluid)
	pp.projectOutOfBnd(flags=gFlags, bnd=params['bnd']+params['dx']*0.5)
	markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	markIsolatedFluidCell(flags=gFlags, mark=FlagEmpty)

	# update velocity of the Lagrangian particles
	updateVelocityFromDeltaPos(parts=pp, vel=pV, x_prev=pVtmp, dt=s.timestep, ptype=pT, exclude=FlagObstacle|FlagFluid)

	# NOTE: We don't need to solve the pressure for isolated cells.
	setPartType(parts=pp, ptype=pT, mark=FlagFluid, stype=FlagEmpty, flags=gFlags, cflag=FlagFluid) # empty -> fluid if they enter again.
	setPartType(parts=pp, ptype=pT, mark=FlagEmpty, stype=FlagFluid, flags=gFlags, cflag=FlagEmpty) # fluid -> empty if they escape

	# keep the valid splashing judgements; the particles may still stay inside the flow (due to a small timestep size)
	curr = np.zeros(pp.pySize(), dtype=dtype_int)
	copyPdataToArrayInt(target=curr, source=pT)
	np_pTimer = np_pTimer - s.timestep
	np_pTimer[(np_pTimer<=0.0)] = 0.0 # time-over of a splashing judgement
	keep = np.zeros(pp.pySize(), dtype=dtype_int)
	keep[(curr==FlagFluid)&(np_pTimer>0.0)] = FlagEmpty # keep valid judgement
	copyArrayToPdataInt(target=pItmp, source=keep)
	pT.setConstIntFlag(s=FlagEmpty, t=pItmp, flag=FlagEmpty) # judgement is still valid -> splashing (empty)

	s.step()

	if output: save_frame(output, s.frame, savingFuncs)
	if 0: gui.screenshot("mlflip_%04d.jpg" % s.frame)
	
tf_sess.close()

stats['fraction'] = float(stats['decision'])/stats['candidate']
stats['reverted'] = stats['decision']-stats['splashed']
if output:
	with open(output+'/run_tf_stats.txt', 'w') as f:
		for key,value in sorted(stats.items()): f.write('{}: {}\n'.format(key, value))

