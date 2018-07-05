/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * GUI extension from python
 *
 ******************************************************************************/

#ifndef _CUSTOMCTRL_H__
#define _CUSTOMCTRL_H__

#include <QSlider>
#include <QLabel>
#include <QCheckBox>
#include <QBoxLayout>
#include "manta.h"

namespace Manta {

// fwd decl.
class Mesh;
class GuiThread;
class MainThread;
	
//! Interface for python declared controls
PYTHON() class CustomControl : public PbClass {
public:
	PYTHON() CustomControl();
	
	virtual void init(QBoxLayout* layout) {};

protected:
};

//! Checkbox with attached text display
class TextCheckbox : public QCheckBox {
	Q_OBJECT
public:
	TextCheckbox(const std::string& name, bool val);
	void attach(QBoxLayout* layout);
	void set(bool v);
	bool get();
	
public slots:
	void update(int v);
		
protected:
	bool mVal;
	QLabel* mLabel;    
	QString mSName;    
};

//! Slider with attached text display
class TextSlider : public QSlider {
	Q_OBJECT
public:
	TextSlider(const std::string& name, float val, float min, float max);
	void attach(QBoxLayout* layout);
	void set(float v);
	float get();
	
public slots:
	void update(int v);
		
protected:
	float mMin, mMax, mScale;
	QLabel* mLabel;    
	QString mSName;    
};
	
//! Links a slider control
PYTHON(name=Slider)
class CustomSlider : public CustomControl  {
public:
	PYTHON() CustomSlider(std::string text, float val, float min, float max);
	virtual void init(QBoxLayout* layout);
	
	PYTHON() float get();
	PYTHON() void set(float v);
	
protected:
	float mMin, mMax, mVal;
	std::string mSName;
	TextSlider* mSlider;
};

//! Links a checkbox control
PYTHON(name=Checkbox)
class CustomCheckbox : public CustomControl  {
public:
	PYTHON() CustomCheckbox(std::string text, bool val);
	virtual void init(QBoxLayout* layout);
	
	PYTHON() bool get();
	PYTHON() void set(bool v);
	
protected:
	bool mVal;
	std::string mSName;
	TextCheckbox* mCheckbox;
};
	

//! GUI adapter class to call from Python
PYTHON() class Gui : public PbClass {
public:
	PYTHON() Gui();
	
	PYTHON() void setBackgroundMesh(Mesh* m);
	PYTHON() void show(bool twoD=false);
	PYTHON() void update();
	PYTHON() void pause();
	PYTHON() PbClass* addControl(PbType t);
	PYTHON() void screenshot(std::string filename);

	// control display upon startup
	PYTHON() void nextRealGrid();
	PYTHON() void nextVec3Grid();
	PYTHON() void nextParts();
	PYTHON() void nextPdata();
	PYTHON() void nextMesh();
	PYTHON() void nextVec3Display();
	PYTHON() void nextMeshDisplay();
	PYTHON() void nextPartDisplay(); 
	PYTHON() void toggleHideGrids();
	PYTHON() void setCamPos(float x, float y, float z);
	PYTHON() void setCamRot(float x, float y, float z);  
	PYTHON() void windowSize(int w, int h);
	
protected:
	GuiThread* mGuiPtr;
	MainThread* mMainPtr;
};
	
} // namespace

#endif

