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

#ifndef _GLWIDGET_H__
#define _GLWIDGET_H__

#include <QGLWidget>
#include <QtOpenGL>
#include "vectorbase.h"

namespace Manta {
	
class GLWidget : public QGLWidget {
Q_OBJECT

public:
	GLWidget(QWidget *parent = NULL);
	~GLWidget();
	
	QSize minimumSizeHint() const;
	QSize sizeHint() const;
	 
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
	void wheelEvent(QWheelEvent *e);     
	void screenshot(QString file);

	void setCamPos(Vec3 pos) { mCamPos = pos; }
	void setCamRot(Vec3 pos) { mRotX = pos.x; mRotY = pos.y; }

public slots:
	void setViewport(const Vec3i& gridsize);
	void keyPressEvent(QKeyEvent* e);
	void keyReleaseEvent(QKeyEvent* e);
	void windowSize(int w, int h);
	 
signals:
	void paintSub();
	void clickLine(QPoint pos, float p0, float p1,float p2, float q0, float q1, float q2);
	void painterEvent(int e, int param=0);
	 
protected:
	bool keyProcess(int key, int mod, bool down);
	void timerEvent(QTimerEvent* e);
	void initializeGL();
	void resizeGL(int w, int h);
	void paintGL();
	void updatePlane(int plane);
	
	enum MoveDir { None = 0, MoveLeft, MoveRight, MoveUp, MoveDown, MoveIn, MoveOut, MoveDirNum };
	
	bool  mMoveState[MoveDirNum];
	bool  mMoveFast;
	QPoint mAnchor, mDownPos;
	Vec3  mCamPos;
	float mRotX, mRotY;
	Vec3i mGridsize;
	int   mPlaneDim, mPlane;
	
	int   mScreenshotNumber;
	int   mWidth, mHeight;
};

} // namespace

#endif
