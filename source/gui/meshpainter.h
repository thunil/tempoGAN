/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Painting mesh objects
 *
 ******************************************************************************/

#ifndef _MESHPAINTER_H_
#define _MESHPAINTER_H_

#include "painter.h"

namespace Manta {
// fwd decl
class Mesh;
	
//! Painter object for Meshes
class MeshPainter : public LockedObjPainter {
	Q_OBJECT
public:
	enum DisplayMode { ModeTrans=0, ModeLines, ModePoints, ModeFlatShade, ModeInvisible, Num_DisplayModes };
	enum BackgroundMode { BModeNormal=0, BModeTrans, BModeInvisible, Num_BackgroundModes };
	enum VorticityMode { VModeFull=0, VModeSmoothed, VModeDiff, VModeSmoke, VModeTex, VModeNone, Num_VorticityModes };
	
	MeshPainter(QWidget* par = 0);
	~MeshPainter();
	
	void paint();
	void attachWidget(QLayout* layout);
	
public slots:
	void setBackgroundMesh(Mesh* bgr);
	
protected:
	std::string getID();    
	void update();
	void updateText();
	void processKeyEvent(PainterEvent e, int param);
	void processSpecificKeyEvent(PainterEvent e, int param);
	void setupLights(bool specular);
	
	Real mColorScale;
	DisplayMode mMode;    
	VorticityMode mVorticityMode;
	BackgroundMode mBackgroundMode;    
	Mesh* mLocalMesh, *mBackground;
	QLabel* mInfo;
	bool mHide;    
};
	
} // namespace

#endif