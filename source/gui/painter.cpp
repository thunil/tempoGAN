/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2016 Tobias Pfaff, Nils Thuerey  
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Base class for objects painting into the GL widget
 *
 ******************************************************************************/

#include "painter.h"
#include "simpleimage.h"
#include <QtOpenGL>
#include <sstream>
#include <iomanip>

using namespace std;

namespace Manta {
	
//******************************************************************************
// Base class
	
void LockedObjPainter::doEvent(int e, int param) {
	// try to obtain valid handle
	 if (!mObject)
		 nextObject();
		
	// filter update events
	if (e == UpdateFull) {
		// always update
		if (mObject) {
			mObject->lock();
			update();
			mObject->unlock();
			mRequestUpdate = false;
		}
	} else if (e == UpdateRequest) {
		// update if resource is available, otherwise wait until next step
		mRequestUpdate = true;        
		if (mObject) {
			if (mObject->tryLock()) {
				update();
				mRequestUpdate = false;
				mObject->unlock();
			}
		}
	} else if (e == UpdateStep) {
		// update if requested only
		if (mRequestUpdate) {
			if (mObject) {
				mObject->lock();
				update();
				mObject->unlock();
				mRequestUpdate = false;
			}
		}
	} else {    
		// pass on all other events
		processKeyEvent((PainterEvent)e, param);
	}
}

void LockedObjPainter::nextObject() {
	if (PbClass::getNumInstances() == 0) return;
	
	int oldIndex = mObjIndex;
	for(;;) {
		mObjIndex = (mObjIndex + 1) % PbClass::getNumInstances();
		if (oldIndex == mObjIndex) break;
		
		PbClass* obj = PbClass::getInstance(mObjIndex);
		if (obj->canConvertTo(getID()) && !obj->isHidden()) {
			mObject = obj;
			doEvent(UpdateRequest);
			return;
		}
		if (oldIndex < 0) oldIndex = 0; // prevent endless loop on first run
	} 
}

//******************************************************************************
// Grid painter

template<class T>
GridPainter<T>::GridPainter(FlagGrid** flags, QWidget* par) 
	: LockedObjPainter(par), mMaxVal(0), mDim(0), mPlane(0), mMax(0), mLocalGrid(NULL), 
	  mFlags(flags), mInfo(NULL), mHide(false), mHideLocal(false), mDispMode(), mValScale()
{
	mDim = 2; // Z plane
	mPlane = 0;
	mInfo = new QLabel(); 
}

template<class T>
GridPainter<T>::~GridPainter() {
	if (mLocalGrid)
		delete mLocalGrid;
}

template<class T>
void GridPainter<T>::attachWidget(QLayout* layout) {
	layout->addWidget(mInfo);
}

template<class T>
void GridPainter<T>::update() {
	Grid<T>* src = (Grid<T>*) mObject;
	
	if (!mLocalGrid) {
		mLocalGrid   = new Grid<T>(src->getParent());
		// int grid is base for resolution
		if (src->getType() & GridBase::TypeInt)
			emit setViewport(src->getSize());
	}
	// reallocate if dimensions changed (or solver)
	if ( (mLocalGrid->getSize() != src->getSize()) || (mLocalGrid->getParent() != src->getParent()) ) { 
		delete mLocalGrid;
		mLocalGrid = new Grid<T>(src->getParent());
		// int grid is base for resolution
		if (src->getType() & GridBase::TypeInt)
			emit setViewport(src->getSize());
	}
	
	mLocalGrid->copyFrom( *src , true ); // copy grid data and type marker
	mLocalGrid->setName(src->getName());
	mLocalGrid->setParent(src->getParent());    
	mMaxVal = mLocalGrid->getMaxAbs();
	
	mPlane = clamp(mPlane, 0, mLocalGrid->getSize()[mDim]-1);
	
	updateText();    
}

template<> string GridPainter<int>::getID()  { return "Grid<int>";  }
template<> string GridPainter<Vec3>::getID() { return "Grid<Vec3>"; }
template<> string GridPainter<Real>::getID() { return "Grid<Real>"; }

template<class T>
void GridPainter<T>::processKeyEvent(PainterEvent e, int param)
{
	if (e == EventSetDim) {
		mDim = param;
		if (!mLocalGrid->is3D()) mDim = 2;
	} else if (e == EventSetMax) {
		mMax = param;
	} else if (e == EventSetPlane) {
		mPlane = param;
		if (mObject) {
			if (mMax>0)
				mPlane = mPlane * mLocalGrid->getSize()[mDim] / mMax;
			mPlane = clamp(mPlane, 0, mLocalGrid->getSize()[mDim]-1);
		}
	} else if (e == EventToggleGridDisplay)
		mHide = !mHide;
	else
		processSpecificKeyEvent(e, param);
	
	updateText();
}

// get scale value for current grid from map, or create new
template<class T>
int GridPainter<T>::getDispMode() {
	if (!mObject || !mLocalGrid) return RealDispOff;
	if (mDispMode.find(mObject) == mDispMode.end()) {
		int dm = RealDispStd; // same for vec & real
		// initialize exceptions, eg levelset
		if (mLocalGrid->getType() & GridBase::TypeLevelset) dm = RealDispLevelset;
		mDispMode[mObject] = dm; 
		return dm;
	} 
	return mDispMode[mObject];
}

template<class T>
void GridPainter<T>::setDispMode(int dm) {
	mDispMode[mObject] = dm;
}

template<class T>
Real GridPainter<T>::getScale() {
	if (!mObject) return 0;
	const int dm = getDispMode();
	std::pair<void*, int> id; id.first=mObject; id.second=dm;
	
	if (mValScale.find(id) == mValScale.end()) {
		// init new scale value
		Real s = 1.0 - VECTOR_EPSILON;
		if (mLocalGrid->getType() & GridBase::TypeVec3)
			s = 0.5 - VECTOR_EPSILON;
		else if (mLocalGrid->getType() & GridBase::TypeLevelset)
			s = 1.0; 
		else if (mLocalGrid->getType() & GridBase::TypeReal) {
			if(dm == RealDispShadeVol ) s = 4.0; // depends a bit on grid size in practice...
			if(dm == RealDispShadeSurf) s = 1.0; 
		}
		mValScale[id] = s;
	}
	return mValScale[id];
	
}

template<class T>
void GridPainter<T>::setScale(Real v) {
	std::pair<void*, int> id; id.first=mObject; id.second=getDispMode();
	mValScale[id] = v;
}

//******************************************************************************
// Grid painter class specializations

template<>
void GridPainter<int>::processSpecificKeyEvent(PainterEvent e, int param) {
	if (e == EventNextInt)
		nextObject();
}

template<>
void GridPainter<Real>::processSpecificKeyEvent(PainterEvent e, int param) {
	if (e == EventNextReal) {
		nextObject();
		mHideLocal = (getDispMode()==VecDispOff); 
	} else if (e == EventScaleRealDown && mObject) {
		setScale( getScale() * 0.5 );
	} else if (e == EventScaleRealUp && mObject) {
		setScale( getScale() * 2.0 );
	} else if (e == EventScaleRealDownSm && mObject) {
		setScale( getScale() * 0.9 );
	} else if (e == EventScaleRealUpSm && mObject) {
		setScale( getScale() * 1.1 );
	} else if (e == EventNextRealDisplayMode) {
		setDispMode( (getDispMode()+1)%NumRealDispModes );
		mHideLocal = (getDispMode()==RealDispOff); 
	}
}

template<>
void GridPainter<Vec3>::processSpecificKeyEvent(PainterEvent e, int param) {
	if (e == EventNextVec) {
		nextObject();
		mHideLocal = (getDispMode()==VecDispOff); 
	} else if (e == EventScaleVecDown && mObject) {
		setScale( getScale() * 0.5 );
	} else if (e == EventScaleVecUp && mObject) {
		setScale( getScale() * 2.0 );
	} else if (e == EventScaleVecDownSm && mObject) {
		setScale( getScale() * 0.9 );
	} else if (e == EventScaleVecUpSm && mObject) {
		setScale( getScale() * 1.1 );
	} else if (e == EventNextVecDisplayMode) {
		setDispMode( (getDispMode()+1)%NumVecDispModes );
		mHideLocal = (getDispMode()==VecDispOff); 
	}
}

template<> void GridPainter<int>::updateText() {
	stringstream s;
	if (mObject) {
		if(mHide) s <<"(hidden) ";
		s << "Int Grid '" << mLocalGrid->getName() << "'" << endl;
	}    
	mInfo->setText(s.str().c_str());    
}

template<> void GridPainter<Real>::updateText() {
	stringstream s;
	
	s << "Display Plane " << mPlane << " [" << (char)('X' + mDim) << "]" << endl << endl;
	if (mObject) {
		s << "Solver '" << mLocalGrid->getParent()->getName() << "'" << endl;
		s << "Grid resolution [" << mLocalGrid->getSizeX() << ", " << mLocalGrid->getSizeY() << ", " << mLocalGrid->getSizeZ() << "]" << endl;
		s << endl;
	}    

	if (mObject && !mHide && !mHideLocal) {
		s << "Real Grid '" << mLocalGrid->getName() << "'" << endl;
		s << "-> Max " << fixed << setprecision(2) << mMaxVal << endl<<"-> Scale " << getScale() << endl;
	}
	mInfo->setText(s.str().c_str());    
}

template<> void GridPainter<Vec3>::updateText() {
	stringstream s;
	if (mObject) {
		s << "Vec Grid '" << mLocalGrid->getName() << "'" << endl;
		if(mHide || mHideLocal) 
			s <<"(hidden) "<< endl;
		else
			s << "-> Max norm " << fixed << setprecision(2) << mMaxVal << endl<<"-> Scale " << getScale() << endl;
	}
	mInfo->setText(s.str().c_str());
}

// compute line intersection with the display plane
Vec3i getQuad(const Vec3& l0, const Vec3& l1, int dim, int plane, Real dx) {
	Vec3 n(0.); n[dim] = 1;
	Vec3 p0 = n*(plane+0.5);
	Vec3 e = (l1-l0)/dx;
	Vec3 e0 = l0/dx;
	Real dotP = dot(p0-e0,n);
	Real dotE = dot(e,n);
	if (dotE == 0) 
		return Vec3i(-1,-1,-1);
	Vec3 s = e0 + (dotP/dotE)*e;
	return toVec3i(s);
}

template<> string GridPainter<int>::clickLine(const Vec3& p0, const Vec3& p1) { 
	if (!mObject) return "";
	Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
	if (!mLocalGrid->isInBounds(s)) return "";
	stringstream m;
	m << "Grid [ " << s.x << ", " << s.y << ", " << s.z << " ]" << endl << mLocalGrid->getName() << ": " << mLocalGrid->get(s) << endl;
	return m.str();
}

template<> string GridPainter<Real>::clickLine(const Vec3& p0, const Vec3& p1) { 
	if (!mObject) return "";
	Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
	if (!mLocalGrid->isInBounds(s)) return "";
	stringstream m;
	m << mLocalGrid->getName() << ": " << setprecision(2) << mLocalGrid->get(s) << endl;
	return m.str();
}

template<> string GridPainter<Vec3>::clickLine(const Vec3& p0, const Vec3& p1) {
	if (!mObject) return "";
	Vec3i s = getQuad(p0,p1,mDim,mPlane,mLocalGrid->getDx());
	if (!mLocalGrid->isInBounds(s)) return "";
	stringstream m;
	m << mLocalGrid->getName() << ": [ " << setprecision(2) << mLocalGrid->get(s).x << ", " <<
										 mLocalGrid->get(s).y << ", " << mLocalGrid->get(s).z << " ]" << endl;
	return m.str();
}


//******************************************************************************
// Actual painting functions

// GL helper functions

// Macro to iterate through one plane
#define FOR_P_SLICE(__g,__dim,__plane) \
	for(Vec3i __g0(__fRange(Vec3i(0,0,0),__dim,__plane)), __g1(__fRange((__g)->getSize(),__dim,__plane+1)), p(__g0); p.z<__g1.z; p.z++) \
		for(p.y=__g0.y; p.y < __g1.y; p.y++) \
			for(p.x=__g0.x; p.x < __g1.x; p.x++)
inline Vec3i __fRange(Vec3i size, int dim, int plane) { Vec3i p(size); p[dim]=plane; return p; }
			  
// coordinate system :
// cell center(i,j,k) -> (i+0.5,j+0.5,k+0.5) / N
// 

void getCellCoordinates(const Vec3i& pos, Vec3 box[4], int dim, bool offset=false) {
	int dim2=(dim+1)%3;
	Vec3 p0(pos.x, pos.y, pos.z);
	Vec3 p1(pos.x+1, pos.y+1, pos.z+1);
	p1[dim] = p0[dim] = pos[dim] + 0.5;

	// display lines with slight offsets
	if(offset) {
		p0 += Vec3(0.01);
		p1 -= Vec3(0.01); 
	}

	box[0] = p0;
	box[3] = p0; box[3][dim2] = p1[dim2];
	box[1] = p1; box[1][dim2] = p0[dim2];
	box[2] = p1;
}
static inline void glVertex(const Vec3& v, const float dx) {
	glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}
void glBox(const Vec3& p0, const Vec3& p1, const float dx) {
	const int box[24] = {0,1,0,2,0,4,7,6,7,5,7,3,1,3,1,5,2,3,2,6,4,5,4,6};
	for (int i=0;i<24;i++) {
		const int b = box[i];
		glVertex(Vec3( (b&1) ? p1.x : p0.x, (b&2) ? p1.y : p0.y, (b&4) ? p1.z : p0.z), dx);
	}
}

// Paint gridlines
template<> void GridPainter<int>::paint() {
	 if (!mObject || mHide || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim])
		return;
	float dx = mLocalGrid->getDx();
	Vec3 box[4];
	glColor3f(0.5,0,0);
	
	bool rbox = true;
	bool skipFluid = mLocalGrid->getSize().max() >= 64; 
	bool drawLines = mLocalGrid->getSize().max() <= 80; 
	if (drawLines) {
		//glDepthFunc(GL_LESS);
		glBegin(GL_LINES);
		FOR_P_SLICE(mLocalGrid, mDim, mPlane) {

			int flag = 0;
			flag = mLocalGrid->get(p);

			if (flag & FlagGrid::TypeObstacle) {
				glColor3f(0.2,0.2,0.2); // dark gray
			} else if (flag & FlagGrid::TypeOutflow) {
				glColor3f(0.9,0.3,0);   // orange
			} else if (flag & FlagGrid::TypeEmpty) {
				glColor3f(0.25,0,0.2);  // dark purple
			} else if (flag & FlagGrid::TypeFluid) {
				if(skipFluid) continue;
				glColor3f(0,0,0.75);    // blue
			} else {
				glColor3f(0.5,0,0); // unknown , medium red
			}

			getCellCoordinates(p, box, mDim, true); 
			for (int n=1;n<=8;n++)
				glVertex(box[(n/2)%4], dx);
		}
		glEnd();
		//glDepthFunc(GL_ALWAYS);        
	}
	
	if (rbox) {
		Vec3 p0(0.0), p1(toVec3(mLocalGrid->getSize())),p(p0);
		glDepthFunc(GL_LESS);
		glBegin(GL_LINES);
		glBox(p0,p1,dx);
		glEnd();
		glDepthFunc(GL_ALWAYS);        
	}
}

// from simpleimage.cpp
void projectImg( SimpleImage& img, const Grid<Real>& val, int shadeMode=0, Real scale=1.);

// Paint box colors
template<> void GridPainter<Real>::paint() {
	if (!mObject || mHide || mHideLocal || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim] || !mFlags || !(*mFlags))
		return;
	
	const int dm     = getDispMode();
	const Real scale = getScale();
	const float dx   = mLocalGrid->getDx();
	Vec3 box[4];
	glBegin(GL_QUADS);

	// "new" drawing style 
	// ignore flags, its a bit dangerous to skip outside info
	if( (dm==RealDispStd) || (dm==RealDispLevelset) ) {

		FOR_P_SLICE(mLocalGrid, mDim, mPlane) 
		{ 
			Real v = mLocalGrid->get(p) * scale; 
			if (dm==RealDispLevelset) {
				v = max(min(v*0.2, 1.0),-1.0);
				if (v>=0)
					glColor3f(v,0,0.5);
				else
					glColor3f(0.5, 1.0+v, 0.);
			} else { // RealDispStd
				if (v>0)
					glColor3f(v,v,v);
				else
					glColor3f(-v,0,0);
			}

			getCellCoordinates(p, box, mDim);
			for (int n=0;n<4;n++) 
				glVertex(box[n], dx);
		}

	}

	if( (dm==RealDispShadeVol) || (dm==RealDispShadeSurf) ) {
		SimpleImage img;

		// note - slightly wasteful, projects all 3 axes!
		int mode = 0;
		if (dm==RealDispShadeSurf) mode = 1;
		projectImg( img, *mLocalGrid, mode, scale );

		FOR_P_SLICE(mLocalGrid, mDim, mPlane) 
		{ 
			Vec3 col(0.); Vec3i s  = mLocalGrid->getSize();

			// "un-transform" projected image
		   	if(mDim==2) col = img.get( 0       + p[0], p[1] );
		   	if(mDim==0) col = img.get( s[0]    + p[2], p[1] );
		   	if(mDim==1) col = img.get( s[0]+s[2]+p[0], p[2] );

			glColor3f(col.x,col.y,col.z); 
			getCellCoordinates(p, box, mDim);
			for (int n=0;n<4;n++) 
				glVertex(box[n], dx);
		}
	}

	glEnd();    
}

// Paint velocity vectors
template<> void GridPainter<Vec3>::paint() {
	if (!mObject || mHide || mHideLocal || mPlane <0 || mPlane >= mLocalGrid->getSize()[mDim])
		return;
	
	const int dm     = getDispMode();
	const Real scale = getScale();
	const float dx   = mLocalGrid->getDx();
	const bool mac   = mLocalGrid->getType() & GridBase::TypeMAC;

	if( (dm==VecDispCentered) || (dm==VecDispStaggered) ) {

		// regular velocity drawing mode
		glBegin(GL_LINES);
			
		FOR_P_SLICE(mLocalGrid, mDim, mPlane) {        
			Vec3 vel = mLocalGrid->get(p) * scale;
			Vec3 pos (p.x+0.5, p.y+0.5, p.z+0.5);
			if (dm==VecDispCentered) {
				if (mac) {
					if (p.x < mLocalGrid->getSizeX()-1) 
						vel.x = 0.5 * (vel.x + scale * mLocalGrid->get(p.x+1,p.y,p.z).x);
					if (p.y < mLocalGrid->getSizeY()-1) 
						vel.y = 0.5 * (vel.y + scale * mLocalGrid->get(p.x,p.y+1,p.z).y);
					if (p.z < mLocalGrid->getSizeZ()-1) 
						vel.z = 0.5 * (vel.z + scale * mLocalGrid->get(p.x,p.y,p.z+1).z);
				}
				glColor3f(0,1,0);
				glVertex(pos, dx);
				glColor3f(1,1,0);
				glVertex(pos+vel*1.2, dx);
			} else if (dm==VecDispStaggered) {
				for (int d=0; d<3; d++) {
					if (fabs(vel[d]) < 1e-2) continue;
					Vec3 p1(pos);
					if (mac)
						p1[d] -= 0.5f;
					Vec3 color(0.0);
					color[d] = 1;
					glColor3f(color.x, color.y, color.z);
					glVertex(p1, dx);
					glColor3f(1,1,0);
					p1[d] += vel[d];
					glVertex(p1, dx);
				}
			}
		}
		glEnd();    
	
	} else if (dm==VecDispUv) {
		// draw as "uv" coordinates (ie rgb), note - this will completely hide the real grid display!
		Vec3 box[4];
		glBegin(GL_QUADS); 
		FOR_P_SLICE(mLocalGrid, mDim, mPlane) 
		{ 
			Vec3 v = mLocalGrid->get(p) * scale; 
			if (mac) {
				if (p.x < mLocalGrid->getSizeX()-1) v.x = 0.5 * (v.x + scale * mLocalGrid->get(p.x+1,p.y,p.z).x);
				if (p.y < mLocalGrid->getSizeY()-1) v.y = 0.5 * (v.y + scale * mLocalGrid->get(p.x,p.y+1,p.z).y);
				if (p.z < mLocalGrid->getSizeZ()-1) v.z = 0.5 * (v.z + scale * mLocalGrid->get(p.x,p.y,p.z+1).z);
			}
			for(int c=0; c<3; ++c) {
				if(v[c]<0.) v[c] *= -1.;
				v[c] = fmod( (Real)v[c], (Real)1.);
			} 
			//v *= mLocalGrid->get(0)[0]; // debug, show uv grid weight as brightness of values
			glColor3f(v[0],v[1],v[2]); 
			getCellCoordinates(p, box, mDim);
			for (int n=0;n<4;n++) 
				glVertex(box[n], dx);
		}
		glEnd();    
	}
}


// explicit instantiation
template class GridPainter<int>;
template class GridPainter<Real>;
template class GridPainter<Vec3>;
	
} // namespace
