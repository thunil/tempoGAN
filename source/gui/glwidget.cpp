/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * QT OpenGL widget
 *
 ******************************************************************************/

#include "glwidget.h"
#ifdef __APPLE__
#   include <OpenGL/glu.h>
#else
#   include <GL/glu.h>
#endif
#include <cmath>
#include "painter.h"

namespace Manta {

GLWidget::GLWidget(QWidget* p): QGLWidget(QGLFormat(QGL::SampleBuffers), p), mRotX(0), mRotY(0), mGridsize(0), mScreenshotNumber(0),
		mWidth(800), mHeight(600)
{
	mPlaneDim = 2; // Y plane
	mPlane = -1;
	mCamPos = Vec3(0, 0, -2);
	for (int i=0; i<MoveDirNum; i++) 
		mMoveState[i] = false;
	mMoveFast = false;
	
	setAutoBufferSwap(true);
	setFocusPolicy(Qt::ClickFocus);
	startTimer(10);
}

GLWidget::~GLWidget()
{

}

QSize GLWidget::minimumSizeHint() const
{
	return QSize(400, 300);
}

QSize GLWidget::sizeHint() const
{
	return QSize(mWidth, mHeight);
}
void GLWidget::windowSize(int w, int h) {
	mWidth = w; mHeight = h;
}

void GLWidget::initializeGL()
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();     
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);    
}

void GLWidget::paintGL()
{
	if (mGridsize.max() == 0) return;
	glDepthFunc(GL_ALWAYS);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(0,0);    
	
	// camera transform
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(mCamPos.x, mCamPos.y, mCamPos.z);
	glRotatef(mRotX,  1.,0.,0.);
	glRotatef(mRotY,  0.,1.,0.);
	Real dx = 1.0 / (Real) mGridsize.max();
	Vec3 sz = toVec3(mGridsize) * (-0.5f * dx);
	
	glTranslatef(sz.x, sz.y, sz.z);
	emit paintSub();
}

void GLWidget::resizeGL(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	GLfloat fov = 45;
	GLfloat zNear = 0.05f;
	GLfloat zFar = 100.0f;
	GLfloat aspect = float(w)/float(h);
	GLfloat fH = tan( float(fov / 360.0f * M_PI) ) * zNear;
	GLfloat fW = fH * aspect;
	glFrustum( -fW, fW, -fH, fH, zNear, zFar );    
	glMatrixMode(GL_MODELVIEW);
	
}

void GLWidget::mouseReleaseEvent(QMouseEvent* event) {
	// only do tooltip if not moving
	QPoint pos = event->pos();
	if ((mDownPos - pos).manhattanLength() == 0) {
		// get GL transform matrices
		int viewport[4];
		GLdouble modelMatrix[16], projMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
		glGetIntegerv(GL_VIEWPORT,viewport);
		
		// obtain click line
		GLdouble line[6], wx=pos.x(), wy=viewport[3]-pos.y();
		if (!gluUnProject(wx,wy,0,modelMatrix,projMatrix,viewport,&line[0],&line[1],&line[2])) return;
		if (!gluUnProject(wx,wy,1.0,modelMatrix,projMatrix,viewport,&line[3],&line[4],&line[5])) return;
		
		// calculate intersection with plane
		emit clickLine(event->globalPos(), line[0],line[1],line[2],line[3],line[4],line[5]);
	}
}

void GLWidget::mouseMoveEvent(QMouseEvent* e)
{
	const float speedRot = 0.2f, speedPan = 0.002f;
	
	QPoint diff = e->pos() - mAnchor;
	 if (e->buttons() & Qt::LeftButton) {
		 mRotX += diff.y() * speedRot;
		 mRotY += diff.x() * speedRot;
		 updateGL();
	 }
	 if (e->buttons() & Qt::RightButton) {
		 mCamPos.x += diff.x() * speedPan;
		 mCamPos.y -= diff.y() * speedPan;
		 updateGL();
	 }

	 mAnchor = e->pos();     
}

void GLWidget::mousePressEvent(QMouseEvent* e)
{
	mDownPos = mAnchor = e->pos();
}

void GLWidget::wheelEvent(QWheelEvent* e)
{
	const float speed = 0.002f;
	mCamPos.z += speed*e->delta();
	updateGL();
}

void GLWidget::timerEvent(QTimerEvent* e)
{
	bool doRepaint = false;

	float speed = 0.005f;
	if (mMoveFast) speed *= 5.;
	
	if (mMoveState[MoveLeft])  { mCamPos.x += speed; doRepaint = true; }
	if (mMoveState[MoveRight]) { mCamPos.x -= speed; doRepaint = true; }
	if (mMoveState[MoveUp])    { mCamPos.y -= speed; doRepaint = true; }
	if (mMoveState[MoveDown])  { mCamPos.y += speed; doRepaint = true; }
	if (mMoveState[MoveOut])   { mCamPos.z -= speed; doRepaint = true; }
	if (mMoveState[MoveIn])    { mCamPos.z += speed; doRepaint = true; }
	if (doRepaint) 
		updateGL();
}

void GLWidget::setViewport(const Vec3i& gridsize) {
	if (mGridsize.x != gridsize.x ||
		mGridsize.y != gridsize.y ||
		mGridsize.z != gridsize.z) {
		if (mPlane < 0) {
			mPlane = gridsize[mPlaneDim] / 2;
		} else {
			Real fac = (Real)gridsize[mPlaneDim] / (Real)mGridsize[mPlaneDim];
			mPlane = (int)(fac * mPlane);
		}
		mGridsize = gridsize;
		emit painterEvent(Painter::EventSetMax, mGridsize[mPlaneDim]);
		emit painterEvent(Painter::EventSetPlane, mPlane);
	}
}

void GLWidget::keyPressEvent(QKeyEvent* e)
{
	if(!keyProcess(e->key(), e->modifiers(), true)) 
		QGLWidget::keyPressEvent(e);
	else 
		updateGL();
}

void GLWidget::keyReleaseEvent(QKeyEvent* e)
{
	if(!keyProcess(e->key(), e->modifiers(), false))
		QGLWidget::keyReleaseEvent(e);
	else
		updateGL();
}

bool GLWidget::keyProcess(int key, int modifier, bool down) 
{
	bool shift = (modifier & Qt::ShiftModifier);
	bool alt   = (modifier & Qt::AltModifier); 
	bool ctrl  = (modifier & Qt::ControlModifier); 
	if      (key == Qt::Key_A) { mMoveState[MoveLeft]  = down; mMoveFast = shift; }
	else if (key == Qt::Key_D) { mMoveState[MoveRight] = down; mMoveFast = shift; }
	else if (key == Qt::Key_W) { mMoveState[MoveUp]    = down; mMoveFast = shift; }
	else if (key == Qt::Key_S) { mMoveState[MoveDown]  = down; mMoveFast = shift; }
	else if (key == Qt::Key_Q) { mMoveState[MoveIn]    = down; mMoveFast = shift; }
	else if (key == Qt::Key_E) { mMoveState[MoveOut]   = down; mMoveFast = shift; }
	else if (down) 
	{
		// only press events
		// note Key_P and Key_L used for play/step in mainwindow.cpp
		if      (key == Qt::Key_Z)                  { /* next "solver" info sometime? */ }
		else if (key == Qt::Key_G)                  { emit painterEvent(Painter::EventToggleGridDisplay); }
		// data grids, first int
		else if (key == Qt::Key_X && shift)         { /* int display mdoes, not yet used */ }
		else if (key == Qt::Key_X)                  { emit painterEvent(Painter::EventNextInt);  updatePlane(mPlane); }
		// real
		else if (key == Qt::Key_C && shift)         { emit painterEvent(Painter::EventNextRealDisplayMode); /* real display modes */ }
		else if (key == Qt::Key_C)                  { emit painterEvent(Painter::EventNextReal); updatePlane(mPlane); } 
		else if (shift && ((key == Qt::Key_Less) ||      
			    (key == Qt::Key_Comma) ) )          { emit painterEvent(Painter::EventScaleRealDownSm); }
		else if (shift && ((key == Qt::Key_Greater) ||
			    (key == Qt::Key_Period) ) )         { emit painterEvent(Painter::EventScaleRealUpSm); }
		else if ((key == Qt::Key_Less) ||      
			    (key == Qt::Key_Comma) )            { emit painterEvent(Painter::EventScaleRealDown); }
		else if ((key == Qt::Key_Greater) ||
			    (key == Qt::Key_Period) )           { emit painterEvent(Painter::EventScaleRealUp); }

		// vec3 grids, scaling can be used with two key combinations (the second one is for international keyboards)
		else if (key == Qt::Key_V && shift)         { emit painterEvent(Painter::EventNextVecDisplayMode); }
		else if (key == Qt::Key_V)                  { emit painterEvent(Painter::EventNextVec);  updatePlane(mPlane); }
		// grid scaling
		else if (key == Qt::Key_Colon)              { emit painterEvent(Painter::EventScaleVecDownSm); }
		else if (key == Qt::Key_QuoteDbl)           { emit painterEvent(Painter::EventScaleVecUpSm); }
		else if (key == Qt::Key_Semicolon)          { emit painterEvent(Painter::EventScaleVecDown); }
		else if (key == Qt::Key_Apostrophe)         { emit painterEvent(Painter::EventScaleVecUp); }

		// particles
		else if (key == Qt::Key_B && shift)         { emit painterEvent(Painter::EventNextParticleDisplayMode); }
		else if (key == Qt::Key_B && alt)           { emit painterEvent(Painter::EventNextSystem); }
		else if (key == Qt::Key_B)                  { emit painterEvent(Painter::EventToggleParticles); }

		else if((key == Qt::Key_ParenLeft) ||       // a bit ugly, but for some reason parentheses dont work in some cases... fall back with dual assignment
			    (key == Qt::Key_9) )                { emit painterEvent(Painter::EventScalePdataDown); }
		else if((key == Qt::Key_ParenRight) ||
			    (key == Qt::Key_0) )                { emit painterEvent(Painter::EventScalePdataUp);   }

		// mesh display
		else if (key == Qt::Key_M && shift)           emit painterEvent(Painter::EventMeshMode);
		else if (key == Qt::Key_BraceLeft )           { emit painterEvent(Painter::EventScaleMeshDown); }
		else if (key == Qt::Key_BracketLeft)          { emit painterEvent(Painter::EventScaleMeshDown); }
		else if (key == Qt::Key_BraceRight)           { emit painterEvent(Painter::EventScaleMeshUp); }
		else if (key == Qt::Key_BracketRight)         { emit painterEvent(Painter::EventScaleMeshUp); }
		// special mesh display modes
		else if (key == Qt::Key_M && alt)             emit painterEvent(Painter::EventMeshColorMode);
		else if (key == Qt::Key_M && ctrl)            emit painterEvent(Painter::EventToggleBackgroundMesh); 
		else if (key == Qt::Key_M)                    emit painterEvent(Painter::EventNextMesh);
		
		// switch display planes
		else if ( (key == Qt::Key_Asterisk) || (key == Qt::Key_8) ) {
			mPlaneDim = (mPlaneDim+1) % 3;            
			emit painterEvent(Painter::EventSetDim, mPlaneDim);
			emit painterEvent(Painter::EventSetMax, mGridsize[mPlaneDim]);
		} 
		// move plane (+shift is faster)
		else if (shift && ((key == Qt::Key_Plus)  || (key == Qt::Key_Equal)) ) { 
			updatePlane(mPlane + 10);
		} else if (shift && ((key == Qt::Key_Minus) || (key == Qt::Key_Underscore)) ) { 
			updatePlane(mPlane - 10);
		} else if (key == Qt::Key_Plus || key == Qt::Key_Equal) { 
			updatePlane(mPlane + 1);
		} else if (key == Qt::Key_Minus) { 
			updatePlane(mPlane - 1);
		}
		else if ( key == Qt::Key_K) {
			QString filename = QString("scr_%1.png").arg(QString::number(mScreenshotNumber), 3, QChar('0'));
			screenshot(filename);
			mScreenshotNumber++;
		}
		
		else return false;
	}
	else return false;
	return true;
}

void GLWidget::screenshot(QString file) {
	grabFrameBuffer().save(file);
}

void GLWidget::updatePlane(int plane) {
	mPlane = clamp(plane, 0, mGridsize[mPlaneDim]);
	emit painterEvent(Painter::EventSetPlane, mPlane);
}



}
