#
# explicit / implicit wave equation solve 
#
import sys
from manta import *
from helperInclude import *

gs  = vec3( 113,127, 1)
s   = Solver(name='main', gridSize = gs, dim=2)

# wave eq settings
implicit   = False
s.timestep = 0.9
cSqr       = 0.12
normalizeMass     = True
useCrankNicholson = False

# allocate grids
h     = s.create(RealGrid)
hprev = s.create(RealGrid)
hnew  = s.create(RealGrid)

curv  = s.create(RealGrid)
vel   = s.create(RealGrid)

flags = s.create(FlagGrid)
flags.initDomain()
flags.fillGrid()

timings = Timings()
if 0 and (GUI):
    gui = Gui()
    gui.show()
    gui.pause()

source = s.create(Box, p0=gs*vec3(0.3,0.3,0.3), p1=gs*vec3(0.5,0.5,0.5))
source.applyToGrid(grid=h,     value=1)
hprev.copyFrom(h)

for t in range(40):

	mass = totalSum( height=h )
	#print "Current mass %f " % mass

	if implicit:
		# implicit solve , cf. 07IntroToPDEs.pdf, page 19
		cgSolveWE( flags=flags, ut=h, utm1=hprev, out=hnew , cSqr=cSqr, crankNic=useCrankNicholson );

	else: 
		# explicit solve , easier-to-read version with explicit velocity integration
		calcSecDeriv2d(h, curv)

		vel.addScaled(curv, cSqr * s.timestep)

		h.addScaled(vel,s.timestep)

		# switch to implicit for second half
		if(t>=20):
			implicit = True


	if normalizeMass:
		normalizeSumTo(h, mass)


	#gui.screenshot( 'out_%04d.png' % t );
	#timings.display()
	s.step()

    
doTestGrid( sys.argv[0], "height" , s, h  , threshold=1e-08 , thresholdStrict=1e-10  )
doTestGrid( sys.argv[0], "vel"    , s, vel, threshold=1e-08 , thresholdStrict=1e-10  )

