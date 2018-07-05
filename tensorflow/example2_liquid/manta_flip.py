# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Fluid implicit particle (FLIP) with surface tension effects (no particle resampling)
#
# ----------------------------------------------------------------------------

import os, argparse, operator, math, random, pickle
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

def normalize(v):
	vmag = math.sqrt(sum(v[i]*v[i] for i in range(len(v))))
	return [ v[i]/vmag for i in range(len(v)) ]

parser = argparse.ArgumentParser(description='FLIP', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o', '--output', default='../data/manta-flip',  help='path to the simulation output')
parser.add_argument('-s', '--seed',   default=1, type=int, help='random seed; used for randomizing initialization')
pargs = parser.parse_known_args()[0]

nogui       = False
pause       = False
output      = os.path.normpath(pargs.output)
savingFuncs = []

# default solver parameters
params                = {}
params['dim']         = 2                  # dimension
params['sres']        = 2                  # sub-resolution per cell
params['dx']          = 1.0/params['sres'] # particle spacing (= 2 x radius)
params['res']         = 64                 # reference resolution
params['len']         = 1.0                # reference length
params['bnd']         = 4                  # boundary cells
params['gref']        = 0                  # real-world gravity
params['jitter']      = 0.0                # jittering particles
params['stref']       = 0.073              # surface tension (reference scale [m]; e.g., 0.073)
params['cgaccuracy']  = 1e-3               # cg solver's threshold
params['fps']         = 30
params['t_end']       = 3.0
params['sdt']         = None
params['frame_saved'] = -1
params['frame_last']  = -1

scaleToManta = float(params['res'])/params['len']
params['gs']    = [params['res']+params['bnd']*2, params['res']+params['bnd']*2, params['res']+params['bnd']*2 if params['dim']==3 else 1]
params['grav']  = params['gref']*scaleToManta
params['stens'] = params['stref']*scaleToManta

random.seed(pargs.seed)

s             = Solver(name='FLIP', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = s.frameLength
s.timestep    = s.frameLength

# prepare grids and particles
gFlags  = s.create(FlagGrid)
_dummyV_ = s.create(MACGrid) # show smaller velocities for UI...
_dummyR_ = s.create(RealGrid) # hide
gV      = s.create(MACGrid)
gVold   = s.create(MACGrid)
gP      = s.create(RealGrid)
gPhi    = s.create(LevelsetGrid)
gCurv   = s.create(RealGrid)

pp      = s.create(BasicParticleSystem)
gIdxSys = s.create(ParticleIndexSystem)
gIdx    = s.create(IntGrid)

pT    = pp.create(PdataInt)
pV    = pp.create(PdataVec3)
pVtmp = pp.create(PdataVec3)

mesh = s.create(name='mesh', type=Mesh) if (params['dim']==3 and not nogui) else None

savingFuncs.append([pp.save, 'particles.uni'])
savingFuncs.append([pV.save, 'particlesVel.uni'])
savingFuncs.append([pT.save, 'particlesType.uni'])

# boundary
gFlags.initDomain(params['bnd']-1)

# fluid
balls, vels = [], []
N_pairs = random.randint(1, 3)
for i in range(N_pairs):
	balls.append([0.5+random.uniform(-0.25, 0), 0.5+random.uniform(-0.25,0.25), 0.5+random.uniform(-0.25,0.25), random.uniform(0.05,0.1)])
	vels.append(list(map(operator.sub, [0.5]*3, balls[-1][:3])))
	balls.append(list(map(operator.add, [0.5]*3, vels[-1])))
	balls[-1].append(random.uniform(0.05,0.1))
	vels.append(list(map(operator.sub, [0.5]*3, balls[-1][:3])))

for i, ball in enumerate(balls):
	ball_c = vec3(ball[0]*params['res']+params['bnd'], ball[1]*params['res']+params['bnd'], ball[2]*params['res']+params['bnd'] if (params['dim']==3) else 0.5)
	obj = s.create(Sphere, center=ball_c, radius=ball[3]*params['res'])
	begin = pp.pySize()
	sampleShapeWithParticles(shape=obj, flags=gFlags, parts=pp, discretization=params['sres'], randomness=params['jitter'], refillEmpty=True, notiming=True)
	end = pp.pySize()
	pT.setConstRange(s=FlagFluid, begin=begin, end=end, notiming=True)
	markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle)

	vel = vels[i]
	if (params['dim']<3): vel[2] = 0
	vel = list(map(operator.mul, normalize(vel), [random.uniform(1, 3)*params['res']]*3))
	pV.setConstRange(s=vec3(vel[0], vel[1], vel[2]), begin=begin, end=end, notiming=True)

	gPhi.join(obj.computeLevelset(), notiming=True)

gui = None
if not nogui:
	gui = Gui()
	gui.show()
	if pause: gui.pause()

if output:
	save_frame(output, s.frame, savingFuncs)
	with open(output+'/params.pickle', 'wb') as f: pickle.dump(params, f)

while (s.timeTotal<params['t_end']): # main loop

	mapPartsToMAC(vel=gV, flags=gFlags, velOld=gVold, parts=pp, partVel=pV, ptype=pT, exclude=FlagEmpty)
	if params['sdt'] is None: s.adaptTimestep(gV.getMax())
	else: s.adaptTimestepByDt(params['sdt'])

	addGravityNoScale(flags=gFlags, vel=gV, gravity=vec3(0, params['grav'], 0))

	gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
	unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)
	getLaplacian(laplacian=gCurv, grid=gPhi)
	gCurv.clamp(-1.0, 1.0)

	if (params['dim']==3 and gui):
		extrapolateLsSimple(phi=gPhi, distance=4, inside=True)
		gPhi.createMesh(mesh)

	setWallBcs(flags=gFlags, vel=gV)
	solvePressure(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'], phi=gPhi, curv=gCurv, surfTens=params['stens'])
	setWallBcs(flags=gFlags, vel=gV)
	extrapolateMACSimple(flags=gFlags, vel=gV)

	# update velocity (general update from FLIP and individual update for Lagrangian particles)
	flipVelocityUpdate(vel=gV, velOld=gVold, flags=gFlags, parts=pp, partVel=pV, flipRatio=0.97, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	addForcePvel(vel=pV, a=vec3(0, params['grav'], 0), dt=s.timestep, ptype=pT, exclude=FlagObstacle|FlagFluid)

	# update position
	pp.getPosPdata(target=pVtmp)
	pp.advectInGrid(flags=gFlags, vel=gV, integrationMode=IntRK4, deleteInObstacle=False, ptype=pT, exclude=FlagObstacle|FlagEmpty)
	eulerStep(parts=pp, vel=pV, ptype=pT, exclude=FlagFluid|FlagObstacle)
	pp.projectOutOfBnd(flags=gFlags, bnd=params['bnd']+params['dx']*0.5)
	markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle)

	# update velocity of the Lagrangian particles
	updateVelocityFromDeltaPos(parts=pp, vel=pV, x_prev=pVtmp, dt=s.timestep, ptype=pT, exclude=FlagFluid|FlagObstacle)

	# NOTE: We don't need to solve the pressure for isolated cells.
	setPartType(parts=pp, ptype=pT, mark=FlagFluid, stype=FlagEmpty, flags=gFlags, cflag=FlagFluid)
	markIsolatedFluidCell(flags=gFlags, mark=FlagEmpty)
	setPartType(parts=pp, ptype=pT, mark=FlagEmpty, stype=FlagFluid, flags=gFlags, cflag=FlagEmpty)

	_dummyV_.copyFrom(gV) 
	_dummyV_.multConst(Vec3(0.02)) 

	s.step()

	if output: save_frame(output, s.frame, savingFuncs)
