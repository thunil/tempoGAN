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

#ifndef _PARTICLEPAINTER_H_
#define _PARTICLEPAINTER_H_

#include "painter.h"
#include "particle.h"

namespace Manta {

//! Painter object for Particle Systems
class ParticlePainter : public LockedObjPainter {
	Q_OBJECT
public:
	ParticlePainter(GridPainter<int>* gridRef, QWidget* par = 0);
	~ParticlePainter();
	
	void paint();
	void attachWidget(QLayout* layout);

	enum PaintModes { PaintOff=0, PaintVel=1, PaintPos=2 };
	
protected:
	std::string getID();
	Real getScale();
	void update();
	void updateText();
	void processKeyEvent(PainterEvent e, int param);
	
	GridPainter<int>* mGridRef;
	ParticleBase* mLocal;
	QLabel* mInfo;

	//! for standard particle systems, this is used to select the data channel (display is changed with displayMode below)
	int mMode; 
	//! for BasicParticleSystem , change the way data channels are displayed
	int mDisplayMode; 

	int mLastPdata;
	bool mHavePdata;
	Real mMaxVal;
	std::string mPdataInfo;
	std::map<PbClass*, Real> mValScale;

	void paintBasicSys();
};    
	
} // namespace

#endif
