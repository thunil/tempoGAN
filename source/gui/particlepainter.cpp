/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Painting particle systems
 *
 ******************************************************************************/

#include <ctime>
#include "particlepainter.h"
#include <sstream>
#include <iomanip>
#include <QtOpenGL>
#include "vortexpart.h"
#include "turbulencepart.h"

using namespace std;

namespace Manta {

ParticlePainter::ParticlePainter(GridPainter<int>* gridRef, QWidget* par) 
	: LockedObjPainter(par), mGridRef(gridRef), mLocal(0), mMode(PaintVel), mDisplayMode(0),
	mLastPdata(-1), mHavePdata(false), mMaxVal(0.)
{    
	mInfo = new QLabel();
}

ParticlePainter::~ParticlePainter() {
	if (mLocal)
		delete mLocal;
}

void ParticlePainter::attachWidget(QLayout* layout) {
	layout->addWidget(mInfo);
}

void ParticlePainter::update() {
	ParticleBase* src = (ParticleBase*) mObject;
	
	// always reallocate
	if (mLocal) 
		delete mLocal;
	
	mLocal = src->clone();
	
	updateText();    
}

string ParticlePainter::getID() { return "ParticleBase"; }

Real ParticlePainter::getScale() {
	if (!mObject) return 0;
	
	if (mValScale.find(mObject) == mValScale.end()) {
		Real s = 1.0;
		//if (mLocalGrid->getType() & GridBase::TypeVec3) s = 0.4;
		mValScale[mObject] = s;
	}
	return mValScale[mObject];
	
}

void ParticlePainter::processKeyEvent(PainterEvent e, int param) {
	if (e == EventNextSystem)
		nextObject();
	else if (e == EventScalePdataDown && mObject)
		mValScale[mObject] = getScale() * 0.5;
	else if (e == EventScalePdataUp && mObject)
		mValScale[mObject] = getScale() * 2.0;
	else if (e == EventToggleParticles) {
		mMode++;  // apply modulo later depending on particle system
		//if(mMode>PaintVel) mMode=PaintOff;
	}
	else if (e == EventNextParticleDisplayMode) {
		mDisplayMode++;
	}
	else return;
		
	updateText();
}

void ParticlePainter::updateText() {
	ostringstream s;
	
	if (mObject && !(mMode==PaintOff) ) {
		s << mLocal->infoString() << endl;
		s << mPdataInfo;
		if(mHavePdata) {
			s << "-> Max " << fixed << setprecision(2) << mMaxVal << "  Scale " << getScale() << endl;
		}
	}
	mInfo->setText( s.str().c_str() );    
}


static inline void glVertex(const Vec3& v, Real dx) {
	glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}

static inline void glColor(const Vec3& color) {
	glColor3f( std::max(0.0f,std::min(1.0f,(float)color.x)), 
			   std::max(0.0f,std::min(1.0f,(float)color.y)), 
			   std::max(0.0f,std::min(1.0f,(float)color.z)) );
}

void ParticlePainter::paint() {
	if (!mObject) return;
	if (mMode == PaintOff) return;
	float dx = mLocal->getParent()->getDx();
	mHavePdata = false;
	mMaxVal = 0.;
	
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST); // disable depth test for particles, clashes with display plane for regular ones
	glDisable(GL_LIGHTING);
	
	// draw points
	if(mLocal->getType() == ParticleBase::VORTEX) {
		VortexParticleSystem* vp = (VortexParticleSystem*) mLocal;
		glColor3f(1,1,0);
		for(int i=0; i<vp->size(); i++) {
			if (vp->isActive(i)) {
				Vec3 pos = (*vp)[i].pos;
			
				glPointSize((*vp)[i].sigma);

				glBegin(GL_POINTS);
				glVertex(pos, dx);
				glEnd();
			}
		}        
	} else if (mLocal->getType() == ParticleBase::FILAMENT) {
		// Filaments don't work yet
		/*VortexFilamentSystem* fp = (VortexFilamentSystem*) mLocal;
		glColor3f(1,1,0);
			
		for(int i=0; i<fp->segSize(); i++) {
			if (!fp->isSegActive(i)) continue;
			const VortexRing& r = fp->seg(i);
			
			glPointSize(1.0);
			glBegin(GL_LINES);
			for(int j=0; j<r.size(); j++) {
				glVertex( (*fp)[r.idx0(j)].pos, dx);
				glVertex( (*fp)[r.idx1(j)].pos, dx);
			}
			glEnd();
		}   */
	} else if(mLocal->getType() == ParticleBase::TURBULENCE) {
		TurbulenceParticleSystem* vp = (TurbulenceParticleSystem*) mLocal;
		glPointSize(2.5);
		glColor3f(0,1,0);
		glBegin(GL_POINTS);
		for(int i=0; i<(int)vp->size(); i++) {
			Vec3 pos = (*vp)[i].pos;
			glColor((*vp)[i].color);
			glVertex(pos, dx);
			
		}   
		glEnd();
		
	} else if(mLocal->getType() == ParticleBase::PARTICLE) {
		paintBasicSys();
	}

	glPointSize(1.0);
	glEnable(GL_DEPTH_TEST); 
}

void ParticlePainter::paintBasicSys() {
	BasicParticleSystem* bp = (BasicParticleSystem*) mLocal;
	//int dim = mGridRef->getDim();
	
	// obtain current plane & draw settings
	int dim = mGridRef->getDim();
	Real factor = mGridRef->getMax() / mLocal->getParent()->getGridSize()[dim];
	int plane = factor * mGridRef->getPlane();
	Real scale = getScale(); 
	float dx = mLocal->getParent()->getDx();

	// draw other particle data, if available
	int pdataId = mMode % (bp->getNumPdata() + 2);
	std::ostringstream infoStr;
	bool drewPoints = false;

	if( pdataId==0 ) {
		// dont draw any points
		infoStr << "Off\n";
		drewPoints = true;
	} else if( pdataId==1 ) {
		// dont draw data, only flags with center below
		infoStr << "Drawing center & flags\n";
	} else if (bp->getNumPdata() > 0)  {
		int pdNum = pdataId-2; // start at 0
		ParticleDataBase* pdb = bp->getPdata(pdNum);

		switch (pdb->getType() ) {

		case ParticleDataBase::TypeReal: {
			ParticleDataImpl<Real>* pdi = dynamic_cast<ParticleDataImpl<Real>*>(pdb);
			if(!pdi) break;
			mHavePdata = true;
			drewPoints = true;
			glPointSize(1.5);
			glBegin(GL_POINTS); 
			for(int i=0; i<(int)bp->size(); i++) {
				if (!bp->isActive(i)) continue;
				Vec3 pos = (*bp)[i].pos; 
				if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
				mMaxVal = std::max( pdi->get(i), mMaxVal );
				Real val = pdi->get(i) * scale;
				glColor3f(0,val,0);
				glVertex(pos, dx); 
			}   
			glEnd();
			infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", real\n";
			} break;

		case ParticleDataBase::TypeInt: {
			ParticleDataImpl<int>* pdi = dynamic_cast<ParticleDataImpl<int>*>(pdb);
			if(!pdi) break;
			mHavePdata = true;
			drewPoints = true;
			glPointSize(1.5);
			glBegin(GL_POINTS); 
			for(int i=0; i<(int)bp->size(); i++) {
				if (!bp->isActive(i)) continue;
				Vec3 pos = (*bp)[i].pos; 
				if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
				Real val = pdi->get(i);
				mMaxVal = std::max( val, mMaxVal );
				val *= scale;
				glColor3f(0,val,0);
				glVertex(pos, dx); 
			}   
			glEnd();
			infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", int\n";
			} break;

		case ParticleDataBase::TypeVec3: { 
			ParticleDataImpl<Vec3>* pdi = dynamic_cast<ParticleDataImpl<Vec3>*>(pdb);
			if(!pdi) break;
			mHavePdata = true;

			// particle vector data can be drawn in different ways...
			mDisplayMode = mDisplayMode%3;

			switch(mDisplayMode) {
			case 0: // lines
				glBegin(GL_LINES); 
				for(int i=0; i<(int)bp->size(); i++) {
					if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
					mMaxVal = std::max( norm(pdi->get(i)), mMaxVal );
					Vec3 val = pdi->get(i) * scale;
					glColor3f(0.5,0.0,0);
					glVertex(pos, dx); 
					pos += val;
					glColor3f(0.5,1.0,0);
					glVertex(pos, dx); 
				}   
				glEnd();
				break;
			case 1:
				// colored points
				glPointSize(2.0);
				glBegin(GL_POINTS); 
				for(int i=0; i<(int)bp->size(); i++) {
					if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
					mMaxVal = std::max( norm(pdi->get(i)), mMaxVal );
					Vec3 val = pdi->get(i) * scale;
					for(int c=0; c<3; ++c) val[c] = fmod( (Real)val[c], (Real)1.);

					glColor3f(val[0],val[1],val[2]);
					glVertex(pos, dx); 
					//pos += val;
					//glColor3f(0.5,1.0,0);
					//glVertex(pos, dx); 
				}   
				glEnd();
				drewPoints = true;
				break;
			case 2:
				glClear(GL_DEPTH_BUFFER_BIT);
				glEnable(GL_DEPTH_TEST); 

				// colored by magnitude all
				glPointSize(2.0);
				glBegin(GL_POINTS); 
				for(int i=0; i<(int)bp->size(); i++) {
					if (!bp->isActive(i)) continue;
					Vec3 pos = (*bp)[i].pos; 
					mMaxVal = std::max( norm(pdi->get(i)), mMaxVal );
					Vec3 val = Vec3( norm( pdi->get(i) * scale ) );
					val[2] += 0.5; // base blue
					for(int c=0; c<3; ++c) val[c] = std::min( (Real)val[c], (Real)1.);

					glColor3f(val[0],val[1],val[2]);
					glVertex(pos, dx); 
				}   
				glEnd();
				drewPoints = true;
				break;
			}

			infoStr << "Pdata '"<<pdi->getName()<<"' #"<<pdNum<<", vec3\n";
			} break;

		default: {
				// skip...
			} break;
		}
	}

	mPdataInfo = infoStr.str(); 
	// enforce refresh upon change
	if(mLastPdata!=pdataId) {
		mLastPdata = pdataId;
		updateText();
	}

	// otherwise draw center
	if(!drewPoints) {
		glPointSize(1.5);
		glBegin(GL_POINTS);

		for(int i=0; i<(int)bp->size(); i++) {
			Vec3 pos = (*bp)[i].pos;
			if (pos[dim] < plane || pos[dim] > plane + 1.0f) continue;
			
			if(!bp->isActive(i) ) {
				glColor3f(1.0, 0., 0.); // deleted, red
			} else if(bp->getStatus(i) & ParticleBase::PNEW ) {
				glColor3f(0.0, 1.0, 0.); // new, green
			} else {
				//glColor3f(0, 0.0, 1.0); // regular, blue
				glColor3f(1.0, 1.0, 1.0); // regular, white - hi contrast
			}
			glVertex(pos, dx);
			
		}   
		glEnd();
	}
	
	// draw basic part sys done
}

} // namespace

