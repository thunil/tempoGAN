# ----------------------------------------------------------------------------
#
# MantaFlow fluid solver framework
# Copyright 2017 Kiwon Um, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Training data generator: Read a simulation data (particle/velocity),
# re-simulate using FLIP in a target resolution, and generate training data
#
# ----------------------------------------------------------------------------

from __future__ import print_function # for calling print(..., end='') in python2
assertNumpy()

import os, glob, sys, argparse, pickle
parser = argparse.ArgumentParser(description='Generate Training Data', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(      '--nogui',  action='store_true',           help='no GUI')
parser.add_argument(      '--pause',  action='store_true',           help='pause')
parser.add_argument('-d', '--dscale', default=0.5, type=float,       help='target resolution (usally, for down-scaling)')
parser.add_argument('-w', '--window', default=1, type=int,           help='window size for sampling features; 1 (default) means 3x3, 2 means 5x5, so on.')
parser.add_argument('-p', '--pfile',  default='particles.uni',       help='file name for particle data')
parser.add_argument('-v', '--vfile',  default='particlesVel.uni',    help='file name for particle velocity data')
parser.add_argument('-t', '--tfile',  default='particlesType.uni',   help='file name for particle type data')
parser.add_argument('-i', '--simdir', default='../data/manta-flip',     help='input path for simulation example; will search for simdir/params.pickle')
parser.add_argument('-o', '--output', default='../data/manta-flip/training_data',  help='output path for training data')
parser.add_argument(      '--prefv',  default='fv',                  help='output path prefix for feature vector data (e.g. 00001/fv/)')
pargs = parser.parse_args()
pargs.output = os.path.normpath(pargs.output)
pargs.simdir = os.path.normpath(pargs.simdir)
indirs = sorted(glob.glob(pargs.simdir + '/0*'))

import numpy as np
dtype_real = np.float32         # NOTE: if double precision, use float64
dtype_int  = np.int32           # NOTE: if int in C is 64bits, use int64

onphi  = True
ongeom = False

loadpath = pargs.simdir+'/params.pickle'
if not os.path.isfile(loadpath): sys.exit('Cannot load the example simulation in '+loadpath)
with open(loadpath, 'rb') as f: params = pickle.load(f)
params['frame_saved'] = -1
params['frame_last']  = -1
pause = pargs.pause

# prepare grids and particles
#  accurate simulation (original)
xl_s             = Solver(name='FLIP-GENDATA-ORI', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
xl_s.cfl         = 1
xl_s.frameLength = 1.0/float(params['fps'])
xl_s.timestepMin = 0
xl_s.timestepMax = xl_s.frameLength
xl_s.timestep    = xl_s.frameLength

xl_gFlags = xl_s.create(FlagGrid)
xl_gR     = xl_s.create(IntGrid)

xl_pp = xl_s.create(BasicParticleSystem)
xl_pT = xl_pp.create(PdataInt)

xl_gFlags.initDomain(params['bnd']-1)

#  downscaled simulation (target; no surface tension effect)
params['res']  *= pargs.dscale
params['bnd']  *= pargs.dscale
params['dx']   *= pargs.dscale
params['grav'] *= pargs.dscale
params['gs']    = [params['res']+params['bnd']*2*pargs.dscale, params['res']+params['bnd']*2*pargs.dscale, params['res']+params['bnd']*2*pargs.dscale if params['dim']==3 else 1]
params['stens'] = None
params['stref'] = None

s             = Solver(name='FLIP-GENDATA-TAR', gridSize=vec3(params['gs'][0], params['gs'][1], params['gs'][2]), dim=params['dim'])
s.cfl         = 1
s.frameLength = 1.0/float(params['fps'])
s.timestepMin = 0
s.timestepMax = xl_s.frameLength
s.timestep    = xl_s.frameLength

gFlags  = s.create(FlagGrid)
_dummyV_ = s.create(MACGrid) # show smaller velocities for UI...
_dummyR_ = s.create(RealGrid) # hide
gR      = s.create(IntGrid)
gV      = s.create(MACGrid)
gVold   = s.create(MACGrid)
gP      = s.create(RealGrid)
gPhi    = s.create(LevelsetGrid)
gPhiSld = s.create(LevelsetGrid)
gIdxSys = s.create(ParticleIndexSystem)
gIdx    = s.create(IntGrid)

pp    = s.create(BasicParticleSystem)
pT    = pp.create(PdataInt)
pV    = pp.create(PdataVec3)
pVtmp = pp.create(PdataVec3)

gFlags.initDomain(params['bnd']-1)

sreg      = pow(int(1.0/pargs.dscale), params['dim'])
fv_N_axis = 2*pargs.window + 1
fv_N_stn  = fv_N_axis*fv_N_axis*(fv_N_axis if params['dim']==3 else 1)
fv_N_row  = params['dim']*fv_N_stn + (fv_N_stn if onphi else 0) + (fv_N_stn if ongeom else 0)
fv_vscale = params['len']/float(params['res'])

def save_features(opath):
    fv = np.zeros(pp.pySize()*fv_N_row, dtype=dtype_real)

    off_feature = 0
    extractFeatureVel(fv=fv, N_row=fv_N_row, off_begin=off_feature, p=pp, vel=gV, scale=fv_vscale, ptype=pT, exclude=FlagObstacle, window=pargs.window)
    off_feature = params['dim']*fv_N_stn
    if onphi:
        extractFeaturePhi(fv=fv, N_row=fv_N_row, off_begin=off_feature, p=pp, phi=gPhi, scale=1.0, ptype=pT, exclude=FlagObstacle, window=pargs.window)
        off_feature += fv_N_stn
    if ongeom:
        extractFeatureGeo(fv=fv, N_row=fv_N_row, off_begin=params['dim']*fv_N_stn, p=pp, flag=gFlags, scale=1.0, ptype=pT, exclude=FlagObstacle, window=pargs.window)
        off_feature += fv_N_stn

    fv = fv.reshape((-1, fv_N_row))
    np.savez_compressed(opath, inputs=fv)

def save_new_splashing_particles(o_path, o_t_path, i_t_curr, i_p_curr, i_p_next):
    xl_pp.load(i_p_curr)
    xl_pT.load(i_t_curr)
    markFluidCells(parts=xl_pp, flags=xl_gFlags, ptype=xl_pT, exclude=FlagObstacle)
    setPartType(parts=xl_pp, ptype=xl_pT, mark=FlagFluid, stype=FlagEmpty|FlagFluid, flags=xl_gFlags, cflag=FlagFluid)

    # 1. kill already splashed particles
    getRegionalCounts(r=xl_gR, flags=xl_gFlags, ctype=FlagFluid)
    #markSmallRegions(flags=xl_gFlags, rcnt=xl_gR, mark=FlagEmpty, exclude=FlagObstacle|FlagOpen, th=sreg)
    markSmallRegions(flags=xl_gFlags, rcnt=xl_gR, mark=FlagEmpty, exclude=FlagObstacle, th=sreg)
    setPartType(parts=xl_pp, ptype=xl_pT, mark=0, stype=FlagFluid, flags=xl_gFlags, cflag=FlagEmpty)

    # 2. mark newly splashing particles
    xl_pp.load(i_p_next)
    markFluidCells(parts=xl_pp, flags=xl_gFlags, ptype=xl_pT, exclude=FlagObstacle)
    getRegionalCounts(r=xl_gR, flags=xl_gFlags, ctype=FlagFluid)
    #markSmallRegions(flags=xl_gFlags, rcnt=xl_gR, mark=FlagEmpty, exclude=FlagObstacle|FlagOpen, th=sreg)
    markSmallRegions(flags=xl_gFlags, rcnt=xl_gR, mark=FlagEmpty, exclude=FlagObstacle, th=sreg)
    setPartType(parts=xl_pp, ptype=xl_pT, mark=FlagEmpty, stype=FlagFluid, flags=xl_gFlags, cflag=FlagEmpty)

    # 3. kill meaningless particles (inside the flow body)
    extendRegion(flags=xl_gFlags, region=FlagEmpty, exclude=FlagObstacle, depth=1)
    setPartType(parts=xl_pp, ptype=xl_pT, mark=0, stype=FlagFluid, flags=xl_gFlags, cflag=FlagFluid)

    if o_t_path: xl_pT.save(o_t_path)

    np_arr = np.zeros(xl_pp.pySize(), dtype=dtype_int)
    copyPdataToArrayInt(target=np_arr, source=xl_pT)
    np.savez_compressed(o_path, labels=np_arr.reshape((-1, 1)))

def drop_zdim(data):
    return np.delete(data, -1, 1)

def save_velocity_modification(o_path, i_t_gt, i_p_next, i_p_curr, i_v_next):
    # i_p_next and i_p_curr: ground truth (high-res) scale
    # i_v_next: target (low-res) scale
    xl_pp.load(i_p_next); xl_pp.getPosPdata(pVtmp)
    xl_pp.load(i_p_curr); xl_pp.getPosPdata(pV) # NOTE: pV will be reloaded so it's safe to use; let's save the memory
    pVtmp.sub(pV); pVtmp.multConst(vec3(pargs.dscale)); pVtmp.multConst(vec3(1.0/s.frameLength))
    pV.load(i_v_next); pVtmp.sub(pV) # dv = (x(n+1) - x(n))/dt - v(n+1)
    np_arr = np.zeros(pp.pySize()*3, dtype=dtype_real)
    copyPdataToArrayVec3(target=np_arr, source=pVtmp)
    np_arr = np_arr.reshape((-1, 3)) if params['dim']==3 else drop_zdim(np_arr.reshape((-1, 3)))

    # clear for non-splashing particles
    lb = np.reshape(np.load(i_t_gt)['labels'], (-1))
    np_arr[lb!=FlagEmpty] = 0
    np.savez_compressed(o_path, modvel=np_arr)

gui = None
if not pargs.nogui:
    gui = Gui()
    gui.show()
    if pargs.pause: gui.pause()    

for i, dir_i in enumerate(indirs):
    if i == len(indirs)-1: break
    print('Frame: {}'.format(dir_i))

    dir_c = os.path.normpath(dir_i)
    dir_n = os.path.normpath(indirs[i+1])

    file_p_c, file_p_n = '{}/{}'.format(dir_c, pargs.pfile), '{}/flip.{}'.format(dir_n, pargs.pfile)
    file_v_c, file_v_n = '{}/{}'.format(dir_c, pargs.vfile), '{}/flip.{}'.format(dir_n, pargs.vfile)
    file_t_c, file_t_n = '{}/{}'.format(dir_c, pargs.tfile), '{}/flip.{}'.format(dir_n, pargs.tfile)
    file_p_gt          = '{}/{}'.format(dir_n, pargs.pfile)
    file_fv            = '{}/{}/input.{}'.format(dir_n, pargs.prefv, 'inputs.npz')
    file_lb            = '{}/{}/label.{}'.format(dir_n, pargs.prefv, 'labels.npz')
    file_vm            = '{}/{}/label.{}'.format(dir_n, pargs.prefv, 'modvel.npz')

    dir_fv = '{}/{}'.format(dir_n, pargs.prefv)
    os.path.isdir(dir_fv) or os.makedirs(dir_fv)

    xl_pp.load(file_p_c)
    pp.load(file_p_c)
    pV.load(file_v_c)
    pV.multConst(vec3(pargs.dscale))
    pT.load(file_t_c)
    markFluidCells(parts=pp, flags=gFlags, ptype=pT, exclude=FlagObstacle)

    frame_last = s.frame
    while (frame_last == s.frame): # NOTE: meaningless codes but included for keeping the simulation context

        mapPartsToMAC(vel=gV, flags=gFlags, velOld=gVold, parts=pp, partVel=pV, ptype=pT, exclude=FlagObstacle|FlagEmpty)
        #s.adaptTimestepByDt(s.frameLength) # NOTE: frame-to-frame
        s.adaptTimestep(1.)

        addGravityNoScale(flags=gFlags, vel=gV, gravity=vec3(0, params['grav'], 0))

        gridParticleIndex(parts=pp, flags=gFlags, indexSys=gIdxSys, index=gIdx)
        unionParticleLevelset(parts=pp, indexSys=gIdxSys, flags=gFlags, index=gIdx, phi=gPhi, radiusFactor=1.0)

        setWallBcs(flags=gFlags, vel=gV)
        solvePressure(flags=gFlags, vel=gV, pressure=gP, cgAccuracy=params['cgaccuracy'], phi=gPhi)
        setWallBcs(flags=gFlags, vel=gV)
        extrapolateMACSimple(flags=gFlags, vel=gV)

        # extract data for feature vector at the beginning of each frame
        save_features(opath=file_fv)

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

    pp.save(file_p_n)
    pV.save(file_v_n)
    pT.save(file_t_n)

    # extract label (here, splash vs non-splash)
    file_t_b = '{}/{}'.format(dir_fv, 'particlesType.uni')
    save_new_splashing_particles(o_path=file_lb, o_t_path=file_t_b, i_t_curr=file_t_c, i_p_curr=file_p_c, i_p_next=file_p_gt)
    save_velocity_modification(o_path=file_vm, i_t_gt=file_lb, i_p_next=file_p_gt, i_p_curr=file_p_c, i_v_next=file_v_n)

################################################################################
# generate training data (inputs/labels) for TensorFlow
os.path.isdir(pargs.output) or os.makedirs(pargs.output)
set_0 = {}
set_1 = {}
opath_0 = pargs.output + '/total_p0.npz' # non-splashing
opath_1 = pargs.output + '/total_p1.npz' # splashing
for dir_i in indirs:
    dir_i  = os.path.normpath(dir_i)
    dir_fv = '{}/{}'.format(dir_i, pargs.prefv)

    # inputs/labels
    paths = {
        'inputs': '{}/input.{}'.format(dir_fv, 'inputs.npz'),
        'labels': '{}/label.{}'.format(dir_fv, 'labels.npz'),
        'modvel': '{}/label.{}'.format(dir_fv, 'modvel.npz')
    }

    if not all([os.path.isfile(paths[x]) for x in paths]):
        print('In {}: Incomplete set of files; skipped'.format(dir_i))
        continue

    labels  = np.load(paths['labels'])['labels']
    if np.sum((labels==FlagEmpty).astype(int))==0:
        print('In {}: No splash particle; skipped'.format(dir_i))
        continue

    del_idx = [k for k, v in enumerate(labels) if (v==FlagObstacle) or (v==0)] # {0: unknown, 1: fluid; 2: obstacle; 4: empty/splash}

    inputs = np.load(paths['inputs'])['inputs']
    inputs = np.delete(inputs, del_idx, 0)

    labels = np.delete(labels, del_idx, 0)
    labels = (labels==FlagEmpty).astype(float)
    modvel = np.load(paths['modvel'])['modvel']
    modvel = np.delete(modvel, del_idx, 0)

    p0_idx = [k for k, v in enumerate(labels) if v==0.0]
    p1_idx = [k for k, v in enumerate(labels) if v==1.0]

    set_0['inputs'] = np.concatenate((set_0['inputs'], inputs[p0_idx]), axis=0) if 'inputs' in set_0 else inputs[p0_idx]
    set_0['labels'] = np.concatenate((set_0['labels'], labels[p0_idx]), axis=0) if 'labels' in set_0 else labels[p0_idx]
    set_0['modvel'] = np.concatenate((set_0['modvel'], modvel[p0_idx]), axis=0) if 'modvel' in set_0 else modvel[p0_idx]
    set_1['inputs'] = np.concatenate((set_1['inputs'], inputs[p1_idx]), axis=0) if 'inputs' in set_1 else inputs[p1_idx]
    set_1['labels'] = np.concatenate((set_1['labels'], labels[p1_idx]), axis=0) if 'labels' in set_1 else labels[p1_idx]
    set_1['modvel'] = np.concatenate((set_1['modvel'], modvel[p1_idx]), axis=0) if 'modvel' in set_1 else modvel[p1_idx]

    print('In {}: Generated {} tuples'.format(dir_i, labels.shape[0]))

print('\nWriting to {}: {} ... '.format(opath_0, list(map(np.shape, set_0.values()))), end='')
np.savez_compressed(opath_0, **set_0)
print('Done.')
print(  'Writing to {}: {} ... '.format(opath_1, list(map(np.shape, set_1.values()))), end='')
np.savez_compressed(opath_1, **set_1)
print('Done.')

# stats
range_idx = [['vel', [0, fv_N_stn*params['dim']], params['dim']]]
if onphi:  range_idx.append(['phi',  [range_idx[-1][1][1], range_idx[-1][1][1]+fv_N_stn], 1])
if ongeom: range_idx.append(['geom', [range_idx[-1][1][1], range_idx[-1][1][1]+fv_N_stn], 1])

stats = { 'input_min_0': np.amin(set_0['inputs'], axis=0),
          'input_max_0': np.amax(set_0['inputs'], axis=0),
          'input_avg_0': np.mean(set_0['inputs'], axis=0),
          'input_std_0': np.std(set_0['inputs'], axis=0),
          'input_min_1': np.amin(set_1['inputs'], axis=0),
          'input_max_1': np.amax(set_1['inputs'], axis=0),
          'input_avg_1': np.mean(set_1['inputs'], axis=0),
          'input_std_1': np.std(set_1['inputs'], axis=0) }

for i in range_idx:
    stats['{}_min_0'.format(i[0])] = np.amin(stats['input_min_0'][i[1][0]:i[1][1]].reshape((-1, i[2])), axis=0)
    stats['{}_max_0'.format(i[0])] = np.amax(stats['input_max_0'][i[1][0]:i[1][1]].reshape((-1, i[2])), axis=0)
    stats['{}_min_1'.format(i[0])] = np.amin(stats['input_min_1'][i[1][0]:i[1][1]].reshape((-1, i[2])), axis=0)
    stats['{}_max_1'.format(i[0])] = np.amax(stats['input_max_1'][i[1][0]:i[1][1]].reshape((-1, i[2])), axis=0)

np.set_printoptions(linewidth=np.inf)
with open('{}/data_stats.txt'.format(pargs.output), 'w') as f:
    for key,value in sorted(stats.items()):
        f.write('{}: {}\n'.format(key, value))
