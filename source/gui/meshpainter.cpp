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

#include "meshpainter.h"
#include <QtOpenGL>
#include "mesh.h"
#include "vortexsheet.h"

using namespace std;

namespace Manta {

MeshPainter::MeshPainter(QWidget* par) 
	: LockedObjPainter(par), mColorScale(1.0), mMode(ModeFlatShade), mVorticityMode(VModeNone), 
	mBackgroundMode(BModeNormal), mLocalMesh(0) , mBackground(0), mHide(false)
{    
	mInfo = new QLabel();
}

MeshPainter::~MeshPainter() {
	if (mLocalMesh)
		delete mLocalMesh;
}

void MeshPainter::attachWidget(QLayout* layout) {
	layout->addWidget(mInfo);
}

void MeshPainter::update() {
	Mesh* src = (Mesh*) mObject;
	
	// always reallocate
	if (mLocalMesh) 
		delete mLocalMesh;
	
	mLocalMesh = src->clone();
	
	updateText();    
}

string MeshPainter::getID() { return "Mesh"; }


void MeshPainter::processKeyEvent(PainterEvent e, int param)
{
	if (e == EventNextMesh)
		nextObject();
	else if (e == EventMeshMode)
		mMode = (DisplayMode) ((mMode+1) % (int)Num_DisplayModes);
	else if (e == EventScaleMeshUp)
		mColorScale *= 2.0;
	else if (e == EventScaleMeshDown)
		mColorScale /= 2.0;
	else if (e == EventMeshColorMode)
		mVorticityMode = (VorticityMode) ((mVorticityMode+1)%(int)Num_VorticityModes);
	else if (e == EventToggleBackgroundMesh)
		mBackgroundMode = (BackgroundMode) ((mBackgroundMode+1)%(int)Num_BackgroundModes);
	else return;
		
	updateText();
}

void MeshPainter::updateText() {
	stringstream s;
	
	if ( mObject && !mHide && mLocalMesh ) {
		s << "Mesh '" << mLocalMesh->getName() << "' [" << mLocalMesh->numTris() << " tris]" << endl;
		if (mMode == ModeFlatShade) s << "DisplayMode: Flatshade" << endl;
		if (mMode == ModeInvisible) s << "DisplayMode: Invisible" << endl;
		if (mMode == ModeLines) s << "DisplayMode: Wireframe" << endl;
		if (mMode == ModePoints) s << "DisplayMode: Pointset" << endl;
		if (mMode == ModeTrans) s << "DisplayMode: Semi-Transparent" << endl;
		if (mLocalMesh->getType() == Mesh::TypeVortexSheet) {
			if (mVorticityMode==VModeFull) s << "Vorticity: full" << endl;
			if (mVorticityMode==VModeSmoothed) s << "Vorticity: gauss smooth" << endl;
			if (mVorticityMode==VModeDiff) s << "Vorticity: gauss diff" << endl;
			if (mVorticityMode==VModeSmoke) s << "Smoke Density" << endl;
			if (mVorticityMode==VModeTex) s << "Texture coordinates 0" << endl;
			if (mVorticityMode!=VModeNone) s << "  Scale " << 1.0/mColorScale << "" << endl;
		}
	}
	mInfo->setText(s.str().c_str());    
}

void MeshPainter::setBackgroundMesh(Mesh* bgr) {
	mBackground = bgr;
}

void MeshPainter::setupLights(bool specular) {
	// control colors
	float max = 1.0;
	float dim = 0.5;
	float dims = specular ? dim : 0;
	float maxs = specular ? max : 0;
	float amb = 0.1;
	float shininess = 50.;
	dim = 0.5; max = 0.75; amb = 0.25;

	float ambient0[] = {amb, amb, amb, max};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0); 

	float diffuse0[]  = {max, dim, dim, 1.0};
	float specular0[] = {maxs, dims, dims, 1.0};
	float position0[] = { 5., 5., 5., 1.0f };
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0); 
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular0); 
	glLightfv(GL_LIGHT0, GL_POSITION, position0);

	float diffuse1[]  = {dim, max, dim, 1.0};
	float specular1[] = {dims, maxs, dims, 1.0};
	float position1[] = {  5., -5., -5., 1.0f };
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse1); 
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular1); 
	glLightfv(GL_LIGHT1, GL_POSITION, position1);

	float diffuse2[]  = {dim, dim, max, 1.0};
	float specular2[] = {dims, dims, maxs, 1.0};
	float position2[] = {  0.3,  2., -10., 1.0f };
	glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse2); 
	glLightfv(GL_LIGHT2, GL_SPECULAR, specular2); 
	glLightfv(GL_LIGHT2, GL_POSITION, position2);

	float specReflection[] = { dims, dims, dims, 1.0f };
	if(shininess == 0.) specReflection[0] = specReflection[1] = specReflection[2] = 0.;
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection); 
	glMateriali(GL_FRONT, GL_SHININESS, shininess); 
	glEnable(GL_LIGHTING); 
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
}

static inline void glVertex(const Vec3& v, Real dx) {
	glVertex3f(v.x * dx, v.y * dx, v.z * dx);
}
static inline void glNormal(const Vec3& v) {
	glNormal3f(v.x, v.y, v.z);
}            

void MeshPainter::paint() {
	if (!mObject || mHide || !mLocalMesh) return;

	Real dx = mLocalMesh->getParent()->getDx();
	
	bool triColor = (mMode == ModeFlatShade) && (mLocalMesh->getType() == Mesh::TypeVortexSheet) && (mVorticityMode!=VModeNone);
	bool nodeColor = (mLocalMesh->getType() == Mesh::TypeVortexSheet) && (mVorticityMode==VModeTex);
	
	// setup OpenGL lighting and material
	const float isoAlpha = 0.4;  
	setupLights(false);
	glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
	glEnable(GL_COLOR_MATERIAL);
	glDepthFunc(GL_LESS);        
	//glDisable(GL_CULL_FACE); 
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
			
	// draw background        
	if (mBackground != NULL && mBackgroundMode != BModeInvisible) {
		if (mBackgroundMode != BModeTrans) {
			glDisable(GL_BLEND);
			glEnable(GL_DEPTH_TEST);
			glColor3f(0.3,0.3,0.5);        
		} else {
			glEnable(GL_BLEND);
			glDisable(GL_DEPTH_TEST);
			glColor4f(0.6,0.6,0.8,0.5);        
		}        
		glEnable(GL_CULL_FACE);
		glPolygonOffset(1.0, 0.5);
		glBegin(GL_TRIANGLES);        
				
		for(int tri=0; tri<mBackground->numTris(); tri++) {
			Vec3 normal = mBackground->getFaceNormal(tri);
			for (int c=0; c<3; c++) {
				glNormal(normal);
				glVertex(mBackground->getNode(tri,c), dx);                
			}
		}
		glEnd();        
		glPolygonOffset(1., 1.);
		glDisable(GL_CULL_FACE);        
	}
	
	setupLights(true);
	if (mMode == ModeFlatShade) {
		glDisable(GL_BLEND);
		glEnable(GL_DEPTH_TEST);
	} else {
		glEnable(GL_BLEND);
		glDisable(GL_DEPTH_TEST);
	}
	if (triColor)
		glDisable(GL_LIGHTING);
		
	// draw triangles        
	if (mMode == ModeFlatShade || mMode == ModeTrans) 
	{
		glEnable(GL_CULL_FACE);    
		glPolygonOffset(1.0, 0.5);
		glBegin(GL_TRIANGLES);        
		
		const int numTris = (int)mLocalMesh->numTris();
		for(int tri=0; tri<numTris; tri++) {
			if (!nodeColor && triColor) {
				VortexSheetInfo& info = ((VortexSheetMesh*)mLocalMesh)->sheet(tri);
				Vec3 v = info.vorticity;
				if (mVorticityMode == VModeSmoothed) v = info.vorticitySmoothed;
				if (mVorticityMode == VModeDiff) v -= info.vorticitySmoothed;
				if (mVorticityMode == VModeSmoke) v = info.smokeAmount / 20.0f;
			
				Vec3 ca = v * 20.0 * mColorScale;
				Vec3 color = Vec3(fabs(ca.x),fabs(ca.y),fabs(ca.z));  
				glColor3f(color.x, color.y, color.z);
			} else if (mLocalMesh->isTriangleFixed(tri))
				glColor3f(0,1,0);
			else if (mLocalMesh->tris(tri).flags & Mesh::FfMarked)
				glColor3f(1,0,0);
			else
				glColor4f(0.5,0.7,1.0, isoAlpha); // blue-ish
				
			for (int c=0; c<3; c++) {
				if (nodeColor) {
					Vec3 tc = ((VortexSheetMesh*)mLocalMesh)->tex1(mLocalMesh->tris(tri).c[c]);
					//Vec3 tc2 = ((VortexSheetMesh*)mLocalMesh)->tex2(mLocalMesh->tris(tri).c[c]);
					//Vec3 tc = gAlpha*tc1+(1-gAlpha)*tc2;
					tc = mColorScale * (tc / toVec3(mLocalMesh->getParent()->getGridSize()));
					tc = nmod(tc, Vec3(1,1,1));
					glColor3f(tc.x, tc.y ,tc.z);
				}
				glNormal(mLocalMesh->getFaceNormal(tri));
				glVertex(mLocalMesh->getNode(tri,c), dx);
			}
		}
		glEnd();        
		glPolygonOffset(1., 1.);
		glDisable(GL_CULL_FACE);        
	}

	// Disable light setup
	glDisable(GL_BLEND);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_TEXTURE_2D);
	
	// draw mesh lines
	if(mMode == ModeLines) {
		glColor3d(1.0, 0.9, 0.9);
		glLineWidth(1.0);
		glBegin(GL_LINES);
		const int numTris = (int)mLocalMesh->numTris();
		for(int tri=0; tri<numTris; tri++)
			for (int j=5; j<5+6; j++)
				glVertex( mLocalMesh->getNode(tri,(j/2)%3), dx);
		glEnd();
	}
	
	// draw vertex points
	if(mMode == ModePoints) {
		static const Vec3 colorSpecial (0.3, 0.5, 0.2);
		//static const Vec3 colortable[] = { Vec3(0.5), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1) };
	
		glPointSize(2.0);
		glBegin(GL_POINTS);
		const int numNodes = (int)mLocalMesh->numNodes();
		for(int i=0; i<numNodes; i++) {
			Vec3 color(0.5, 0.5, 0.5);
			if (mLocalMesh->isNodeFixed(i))
				color = Vec3(0,1,0);
			else if (mLocalMesh->nodes(i).flags & Mesh::NfMarked)
				color = Vec3(1,0,0);
			//int flags = mLocalMesh->flags(i);
			
			glColor3f(color.x, color.y, color.z);
			glVertex(mLocalMesh->nodes(i).pos, dx);
		}
		glEnd();
		glPointSize(1.0);
	}    
}

	
	
} // namespace
