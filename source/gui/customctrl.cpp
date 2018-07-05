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

#include "customctrl.h"
#include "qtmain.h"
#include "painter.h"

using namespace std;

namespace Manta {
	
// *****************************************************************************
// Slider class

CustomControl::CustomControl() : PbClass(0) {
}

CustomCheckbox::CustomCheckbox(string name, bool val) : mVal(val), mSName(name), mCheckbox(0) {
}

void CustomCheckbox::init(QBoxLayout* layout) {
	mCheckbox = new TextCheckbox(mSName, mVal);
	QObject::connect(mCheckbox, SIGNAL(stateChanged(int)), mCheckbox, SLOT(update(int)));
	mCheckbox->attach(layout);
}

bool CustomCheckbox::get() {
	if (!mCheckbox) throw Error("Slider is not attached yet!");
	return mCheckbox->get();
}
void CustomCheckbox::set(bool v) {
	if (!mCheckbox) throw Error("Slider is not attached yet!");
	mCheckbox->set(v);
}

CustomSlider::CustomSlider(string name, float val, float min, float max) : 
	mMin(min), mMax(max), mVal(val), mSName(name), mSlider(0)
{
}

void CustomSlider::init(QBoxLayout* layout) {
	mSlider = new TextSlider(mSName, mVal, mMin, mMax);
	QObject::connect(mSlider, SIGNAL(valueChanged(int)), mSlider, SLOT(update(int)));
	mSlider->attach(layout);
}

float CustomSlider::get() {
	if (!mSlider) throw Error("Slider is not attached yet!");
	return mSlider->get();
}

void CustomSlider::set(float v) {
	if (!mSlider) throw Error("Slider is not attached yet!");
	mSlider->set(v);
}

TextSlider::TextSlider(const string& name, float val, float vmin, float vmax) : 
	QSlider(Qt::Horizontal), mMin(vmin), mMax(vmax), mSName(name.c_str())
{
	mLabel = new QLabel();
	mScale = 1000;
	setMinimum(0);
	setMaximum(999);    
	set(val);
	update(0);
 }

void TextSlider::attach(QBoxLayout* layout) {
	layout->addWidget(mLabel);
	layout->addWidget(this);    
}

void TextSlider::update(int val) {
	float v = get();
	QString num;
	num.sprintf("%.2g", v);
	mLabel->setText(mSName + ":  " + num);    
}

float TextSlider::get() {
	float va = mMin + (mMax-mMin) / mScale * (float)value();
	return clamp(va, mMin, mMax);
}

void TextSlider::set(float v) {
	float va = clamp(v, mMin, mMax);
	va = (va - mMin) / (mMax-mMin) * mScale;
	setValue((int)(va+0.5));
}

TextCheckbox::TextCheckbox(const string& name, bool val) : 
	QCheckBox(), mVal(val), mSName(name.c_str())
{
	mLabel = new QLabel();
	set(val);
	mLabel->setText(mSName);
 }

void TextCheckbox::attach(QBoxLayout* layout) {
	QLayout* lay = new QHBoxLayout;    
	lay->setAlignment(Qt::AlignLeft);
	lay->addWidget(this);
	lay->addWidget(mLabel);
	layout->addLayout(lay);
}

void TextCheckbox::update(int val) {
}

bool TextCheckbox::get() {
	return isChecked();
}

void TextCheckbox::set(bool v) {
	setChecked(v);
}


	
// **************************************************************************************
// GUI class

extern MainThread* gMainThread;
extern GuiThread* gGuiThread;

Gui::Gui() : 
	PbClass(NULL), mGuiPtr(gGuiThread), mMainPtr(gMainThread) {     
	if( getenv("MANTA_DISABLE_UI") && atoi( getenv("MANTA_DISABLE_UI") )) { errMsg("Cannot create GUI object, GUI disabled"); }
}

void Gui::setBackgroundMesh(Mesh* m) {
	mGuiPtr->getWindow()->setBackground(m);
}
void Gui::show(bool twoD) {
	if(twoD)
		mMainPtr->send( (int)MainWnd::EventSet2DCam );
	mMainPtr->sendAndWait( (int)MainWnd::EventGuiShow );
}
void Gui::update() { 
	updateQtGui(true, -1,-1., "");
}
void Gui::pause() {
	mMainPtr->sendAndWait((int)MainWnd::EventFullUpdate);         
	mGuiPtr->getWindow()->pause();         
}
void Gui::screenshot(string filename) {
	QString s(filename.c_str());
	QMetaObject::invokeMethod(mGuiPtr->getWindow(), "screenshot", Q_ARG(QString, s));    
}

void Gui::nextRealGrid() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextRealGrid" ); }
void Gui::nextVec3Grid() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextVec3Grid" ); }
void Gui::nextParts()    { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextParts" ); }
void Gui::nextPdata()    { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextPdata" ); }
void Gui::nextMesh()     { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextMesh" ); }

void Gui::nextVec3Display() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextVec3Display" ); }
void Gui::nextPartDisplay() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextPartDisplay" ); }
void Gui::nextMeshDisplay() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "nextMeshDisplay" ); }

void Gui::toggleHideGrids() { QMetaObject::invokeMethod(mGuiPtr->getWindow(), "toggleHideGrids" ); }

void Gui::setCamPos(float x, float y, float z) {
	QMetaObject::invokeMethod(mGuiPtr->getWindow(), "setCamPos", Q_ARG(float, x), Q_ARG(float, y), Q_ARG(float, z));    
}
void Gui::setCamRot(float x, float y, float z) {
	QMetaObject::invokeMethod(mGuiPtr->getWindow(), "setCamRot", Q_ARG(float, x), Q_ARG(float, y), Q_ARG(float, z));    
}
void Gui::windowSize(int w, int h) {
	QMetaObject::invokeMethod(mGuiPtr->getWindow(), "windowSize", Q_ARG(int, w), Q_ARG(int, h));    
}

PbClass* Gui::addControl(PbType t) {
	_args.add("nocheck",true);
	if (t.str() == "")
		throw Error("Need to specify object type. Use e.g. gui.create(Slider, ...)");
	
	PbClass* obj = PbClass::createPyObject(t.str(), "", _args, this);
	if (!obj || !obj->canConvertTo("CustomControl"))
		throw Error("gui.create() can only create CustomControl-based objects");
	
	QMetaObject::invokeMethod(gGuiThread->getWindow(), "addControl", Q_ARG(void*, (void*)obj));    
	
	return obj;
}


} // namespace
