/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011-2016 Tobias Pfaff, Nils Thuerey  
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Loading and writing grids and meshes to disk
 *
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#if NO_ZLIB!=1
extern "C" { 
#include <zlib.h>
}
#endif

#include "mantaio.h"
#include "grid.h"
#include "mesh.h"
#include "vortexsheet.h"
#include  <cstring>

using namespace std;

namespace Manta {


//*****************************************************************************
// mesh data
//*****************************************************************************

void readBobjFile(const string& name, Mesh* mesh, bool append) {
	debMsg( "reading mesh file " << name ,1);
	if (!append)
		mesh->clear();
	else
		errMsg("readBobj: append not yet implemented!");

#	if NO_ZLIB!=1
	const Real dx = mesh->getParent()->getDx();
	const Vec3 gs = toVec3( mesh->getParent()->getGridSize() );

	gzFile gzf = gzopen(name.c_str(), "rb1"); // do some compression
	if (!gzf)
		errMsg("readBobj: unable to open file");
	
	// read vertices
	int num = 0;
	gzread(gzf, &num, sizeof(int));
	mesh->resizeNodes(num);
	debMsg( "read mesh , verts "<<num,1);
	for (int i=0; i<num; i++) {
		Vector3D<float> pos;
		gzread(gzf, &pos.value[0], sizeof(float)*3);
	   	mesh->nodes(i).pos = toVec3(pos);

		// convert to grid space
		mesh->nodes(i).pos /= dx;
		mesh->nodes(i).pos += gs*0.5;
	}
	
	// normals
	num = 0;
	gzread(gzf, &num, sizeof(int));
	for (int i=0; i<num; i++) {
		Vector3D<float> pos;
		gzread(gzf, &pos.value[0], sizeof(float)*3);
	   	mesh->nodes(i).normal = toVec3(pos);
	}
	
	// read tris
	num = 0;
	gzread(gzf, &num, sizeof(int));
	mesh->resizeTris( num );
	for(int t=0; t<num; t++) {
		for(int j=0; j<3; j++) { 
			int trip = 0;
			gzread(gzf, &trip, sizeof(int)); 
			mesh->tris(t).c[j] = trip;
		}
	} 
	// note - vortex sheet info ignored for now... (see writeBobj)
	gzclose( gzf );    
	debMsg( "read mesh , triangles "<<mesh->numTris()<<", vertices "<<mesh->numNodes()<<" ",1 );
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

void writeBobjFile(const string& name, Mesh* mesh) {
	debMsg( "writing mesh file " << name ,1);
#	if NO_ZLIB!=1
	const Real  dx = mesh->getParent()->getDx();
	const Vec3i gs = mesh->getParent()->getGridSize();
	
	gzFile gzf = gzopen(name.c_str(), "wb1"); // do some compression
	if (!gzf)
		errMsg("writeBobj: unable to open file");
	
	// write vertices
	int numVerts = mesh->numNodes();
	gzwrite(gzf, &numVerts, sizeof(int));
	for (int i=0; i<numVerts; i++) {
		Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
		// normalize to unit cube around 0
		pos -= toVec3f(gs)*0.5;
		pos *= dx;
		gzwrite(gzf, &pos.value[0], sizeof(float)*3);
	}
	
	// normals
	mesh->computeVertexNormals();
	gzwrite(gzf, &numVerts, sizeof(int));
	for (int i=0; i<numVerts; i++) {
		Vector3D<float> pos = toVec3f(mesh->nodes(i).normal);
		gzwrite(gzf, &pos.value[0], sizeof(float)*3);
	}
	
	// write tris
	int numTris = mesh->numTris();
	gzwrite(gzf, &numTris, sizeof(int));
	for(int t=0; t<numTris; t++) {
		for(int j=0; j<3; j++) { 
			int trip = mesh->tris(t).c[j];
			gzwrite(gzf, &trip, sizeof(int)); 
		}
	}
	
	// per vertex smoke densities
	if (mesh->getType() == Mesh::TypeVortexSheet) {
		VortexSheetMesh* vmesh = (VortexSheetMesh*) mesh;
		int densId[4] = {0, 'v','d','e'};
		gzwrite(gzf, &densId[0], sizeof(int) * 4); 

		// compute densities
		vector<float> triDensity(numTris);
		for (int tri=0; tri < numTris; tri++) {
			Real area = vmesh->getFaceArea(tri);
			if (area>0)
				triDensity[tri] = vmesh->sheet(tri).smokeAmount;
		}
		
		// project triangle data to vertex
		vector<int> triPerVertex(numVerts);
		vector<float> density(numVerts);
		for (int tri=0; tri < numTris; tri++) {
			for (int c=0; c<3; c++) {
				int vertex = mesh->tris(tri).c[c];
				density[vertex] += triDensity[tri];
				triPerVertex[vertex]++;
			}
		}
		
		// averaged smoke densities
		for(int point=0; point<numVerts; point++) {
			float dens = 0;
			if (triPerVertex[point]>0)
				dens = density[point] / triPerVertex[point];
			gzwrite(gzf, &dens, sizeof(float));             
		}
	}
	
	// vertex flags
	if (mesh->getType() == Mesh::TypeVortexSheet) {
		int Id[4] = {0, 'v','x','f'};
		gzwrite(gzf, &Id[0], sizeof(int) * 4); 

		// averaged smoke densities
		for(int point=0; point<numVerts; point++) {
			float alpha = (mesh->nodes(point).flags & Mesh::NfMarked) ? 1: 0;
			gzwrite(gzf, &alpha, sizeof(float));             
		}
	}

	gzclose( gzf );    
#	else
	debMsg( "file format not supported without zlib" ,1);
#	endif
}

void readObjFile(const std::string& name, Mesh* mesh, bool append) {
	ifstream ifs (name.c_str());
	
	if (!ifs.good())
		errMsg("can't open file '" + name + "'");
	
	if (!append)
		mesh->clear();
	int nodebase = mesh->numNodes();
	while(ifs.good() && !ifs.eof()) {
		string id;
		ifs >> id;
		
		if (id[0] == '#') {
			// comment
			getline(ifs, id);
			continue;
		}
		if (id == "vt") {
			// tex coord, ignore            
		} else if (id == "vn") {
			// normals, ignore            
		} else if (id == "v") {
			// vertex
			Node n;
			ifs >> n.pos.x >> n.pos.y >> n.pos.z;
			mesh->addNode(n);
		} else if (id == "g") {
			// group
			string group;
			ifs >> group;
		} else if (id == "f") {
			// face
			string face;
			Triangle t;
			for (int i=0; i<3; i++) {
				ifs >> face;
				if (face.find('/') != string::npos)
					face = face.substr(0, face.find('/')); // ignore other indices
				int idx = atoi(face.c_str()) - 1;
				if (idx < 0)
					errMsg("invalid face encountered");
				idx += nodebase;
				t.c[i] = idx;
			}
			mesh->addTri(t);
		} else {
			// whatever, ignore
		}
		// kill rest of line
		getline(ifs, id);   
	}
	ifs.close();    
}

// write regular .obj file, in line with bobj.gz output (but only verts & tris for now)
void writeObjFile(const string& name, Mesh* mesh) {
	const Real  dx = mesh->getParent()->getDx();
	const Vec3i gs = mesh->getParent()->getGridSize();

	ofstream ofs(name.c_str());
	if (!ofs.good())
		errMsg("writeObjFile: can't open file " << name);

	ofs << "o MantaMesh\n";
	
	// write vertices
	int numVerts = mesh->numNodes();
	for (int i=0; i<numVerts; i++) {
		Vector3D<float> pos = toVec3f(mesh->nodes(i).pos);
		// normalize to unit cube around 0
		pos -= toVec3f(gs)*0.5;
		pos *= dx;
		ofs << "v "<< pos.value[0] <<" "<< pos.value[1] <<" "<< pos.value[2] <<" "<< "\n";
	}
	
	// no normals for now
	
	// write tris
	int numTris = mesh->numTris();
	for(int t=0; t<numTris; t++) {
		ofs << "f "<< (mesh->tris(t).c[0]+1) <<" "<< (mesh->tris(t).c[1]+1) <<" "<< (mesh->tris(t).c[2]+1) <<" "<< "\n";
	}

	ofs.close();
}

} //namespace
