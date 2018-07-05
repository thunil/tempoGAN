/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Adapted from bobj-loader (C) 2008 Johannes Schmidt
 *
 ******************************************************************************/

// NT - global switch to toggle on/off the use
// of the bobj.gz normals. These seem to cause problems
// sometimes, unlocking would help, but I haven't figured
// out how to do this automatically after each load.
static const bool useNormals = false;

#include "zlib.h"
#include <vector>

#include <maya/MFnStringData.h>
#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVector.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnPlugin.h>

#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMeshData.h>

#include <maya/MDagPath.h>
#include <maya/MFloatArray.h>

#include <maya/MIOStream.h>

MStatus returnStatus;

// debugging manually on/off
#ifndef DODEBUG
#define DODEBUG 1
#endif // DEBUG

#if defined(_DEBUG)||defined(DODEBUG)
#define debMsg(x) std::cerr<<"bobjFluidObject: "<<x<<std::endl;
#else
#define debMsg(x)
#endif


#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}


#define bailOut(x)                                    \
	{ std::cerr << "bobjFluidObject Error: " << x << std::endl; \
	if (filename) delete filename;                    \
	return false; }                                   \

class bobjFluidObject : public MPxNode
{
public:
					bobjFluidObject() {};
	virtual 		~bobjFluidObject() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject time;
	static MObject inFileMask;
	static MObject scale;
	static MObject indexOffset;
	static MObject outputMesh;
	static MTypeId id;

protected:
	bool loadMeshData(const MTime& time, const char* inFileMask, int indexOffset, float scaleFloat, MObject& outData, MStatus& stat);
	void createNullMesh(MObject& outData, MStatus& stat);
};

MObject bobjFluidObject::time;
MObject bobjFluidObject::inFileMask;
MObject bobjFluidObject::indexOffset;
MObject bobjFluidObject::outputMesh;
MObject bobjFluidObject::scale;
MTypeId bobjFluidObject::id( 0x80000 );

void* bobjFluidObject::creator()
{
	return new bobjFluidObject;
}

MStatus bobjFluidObject::initialize()
{
	MFnUnitAttribute unitAttr;
	MFnTypedAttribute typedAttr;
	MFnNumericAttribute numericAttr;
	MFnStringData stringData;
	MObject string;

	MStatus returnStatus;

	/*
	 * Create input attributes
	 */
	bobjFluidObject::time = unitAttr.create( "time", "tm",
										  MFnUnitAttribute::kTime,
										  0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating bobjFluidObject time attribute\n");
	returnStatus = addAttribute(bobjFluidObject::time);
	McheckErr(returnStatus, "ERROR adding time attribute\n");


	string = stringData.create(&returnStatus);
	McheckErr(returnStatus, "ERROR creating string\n");
	bobjFluidObject::inFileMask = typedAttr.create( "inFileMask", "fil",
		MFnData::kString, string,
		&returnStatus);
	McheckErr(returnStatus, "ERROR creating bobjFluidObject inFileMask attribute\n");
	/*MFnStringData inFileMaskData(inFileMask, returnStatus);
	McheckErr(returnStatus, "ERROR retrieving inFileMask data\n");
	inFileMaskData.set(*/
	returnStatus = addAttribute(bobjFluidObject::inFileMask);
	McheckErr(returnStatus, "ERROR adding inFileMask attribute\n");

	bobjFluidObject::indexOffset = numericAttr.create( "indexOffset", "ofs",
		MFnNumericData::kInt,
		-1,
		&returnStatus);
	McheckErr(returnStatus, "ERROR creating bobjFluidObject indexOffset attribute\n");
	returnStatus = addAttribute(bobjFluidObject::indexOffset);
	McheckErr(returnStatus, "ERROR adding indexOffset attribute\n");

	bobjFluidObject::scale = numericAttr.create( "scale", "sc",
		MFnNumericData::kFloat,
		1.0,
		&returnStatus);
	McheckErr(returnStatus, "ERROR creating bobjFluidObject scale attribute\n");
	returnStatus = addAttribute(bobjFluidObject::scale);
	McheckErr(returnStatus, "ERROR adding scale attribute\n");


	/*
	 * Create output attributes
	 */

	bobjFluidObject::outputMesh = typedAttr.create( "outputMesh", "out",
												 MFnData::kMesh,
												 &returnStatus ); 
	McheckErr(returnStatus, "ERROR creating bobjFluidObject output attribute\n");
	typedAttr.setStorable(false);
	returnStatus = addAttribute(bobjFluidObject::outputMesh);
	McheckErr(returnStatus, "ERROR adding outputMesh attribute\n");



	/*
	 * Set dependencies
	 */

	returnStatus = attributeAffects(bobjFluidObject::time,
								    bobjFluidObject::outputMesh);
	McheckErr(returnStatus, "ERROR in making attribute time affect outputMesh\n");
	returnStatus = attributeAffects(bobjFluidObject::scale,
								    bobjFluidObject::outputMesh);
	McheckErr(returnStatus, "ERROR in making attribute scale affect outputMesh\n");
	returnStatus = attributeAffects(bobjFluidObject::inFileMask,
								    bobjFluidObject::outputMesh);
	McheckErr(returnStatus, "ERROR in making attribute inFileMask affect outputMesh\n");

	return MS::kSuccess;
}

void bobjFluidObject::createNullMesh(MObject& outData, MStatus& stat)
{
	MFloatPointArray points;
	MFnMesh meshFS;
	int dummy = 0;

	MObject newMesh = meshFS.create(0, 0, points, 
									dummy, dummy,
									outData, &stat);
}

bool bobjFluidObject::loadMeshData(const MTime& time,
							  const char* inFileMask,
							  int indexOffset,
							  float scale,
							  MObject& outData,
							  MStatus& stat)
{
	int frameNo;
	char* filename = NULL;

	// compute frame index
	frameNo = (int)time.as( MTime::kFilm );
	frameNo += indexOffset;

	// validate input file mask
	if (!inFileMask) 
		bailOut("inFileMask is NULL");
	if (strlen(inFileMask) < 3)
		bailOut("inFileMask is < 3")

	debMsg("bobjFluidObject::loadMeshData: Starting for " << inFileMask<<", frame "<<frameNo );

	const char* perc = strrchr(inFileMask, '%');
	if (perc == NULL)
		bailOut("invalid file mask string " << inFileMask);
	if (!((perc[1]=='d')||(perc[2]=='d')||(perc[3]=='d'))) 
		bailOut("invalid file mask string " << inFileMask);

	// get rid of "pre_" in string? if yes, turn "...pre_..." into "..._..." to load final bobj's
	bool doLoadFinal = false;
	if(getenv("DDF_LOADFINAL")) {
		if(atoi(getenv("DDF_LOADFINAL"))>=1) {
			doLoadFinal=true;
		}
	}

	// get input file name
	filename = new char[strlen(inFileMask)+20];
	if(!doLoadFinal) {
		// normal load
		sprintf(filename, inFileMask, frameNo);
	} else {
		// load final ones
		char *tmpname = new char[strlen(inFileMask)+20];
		strcpy(tmpname, inFileMask);

		char* preloc = strstr(tmpname, "pre_");
		if(preloc) {
			// remove "pre" by copying original end over
			//printf("bobjFluidObject  %s - %s - %s \n", tmpname, preloc, &inFileMask[preloc-tmpname+3]);
			strcpy(preloc, &inFileMask[preloc-tmpname+3]);
			//printf("bobjFluidObject  after %s  \n", tmpname);
		} else {
			debMsg("bobjFluidObject::loadMeshData: error - DDF_LOADFINAL set, but no pre_ in " << filename);
			return false;
		}
		sprintf(filename, tmpname, frameNo);
	}

	debMsg("bobjFluidObject::loadMeshData: Processing file " << filename);

	gzFile gzf;
	int numVerts, numNorms, numTris;
	MFloatPointArray points;
	MFloatVectorArray normals;
	MIntArray faceCounts, faceConnects;
	MFnMesh meshFn;

	// open file
	gzf = gzopen(filename, "rb1");
	if (!gzf)
		bailOut("cannot open file " << filename);
	
	// read vertices
	if (gzread(gzf, &numVerts, 4) != 4)
		bailOut("read error in file " << filename);

	debMsg("Reading " << numVerts << " vertices");
	for (int i = 0; i<numVerts; i++) {
		float v[3];
		for (int j=0; j<3; j++) {
			if (gzread(gzf, &v[j], sizeof(float)) != sizeof(float))
				bailOut("read error in file " << filename);
			v[j] *= scale;
			//if(i<20) debMsg("ReadV "<<i<<","<<j<<" "<<v[j] );
		}
		points.append(MFloatPoint(v[0], v[1], v[2]));
	}

	float pmin[3], pmax[3];
	for(int j=0; j<3; j++) {
		pmin[j] = pmax[j] = 0.;
	}
	if(numVerts>0) {
		for(int j=0; j<3; j++) {
			pmin[j] = points[0][j];
			pmax[j] = points[0][j];
		}
		for(int i=1; i<numVerts; i++) {
			for(int j=0; j<3; j++) {
				if(points[i][j] < pmin[j]) pmin[j] = points[i][j];
				if(points[i][j] > pmax[j]) pmax[j] = points[i][j];
			}
		}
	}
	debMsg("Vertex bounds min=" << pmin[0]<<","<<pmin[1]<<","<<pmin[2]<<"," << " to max="<< pmax[0]<<","<<pmax[1]<<","<<pmax[2]);

	// read normals
	if (gzread(gzf, &numNorms, 4) != 4)
		bailOut("read error in file " << filename);

	debMsg("Reading " << numNorms << " normals");
	for (int i = 0; i<numNorms; i++) {
		float v[3];
		for (int j=0; j<3; j++) {
			if (gzread(gzf, &v[j], sizeof(float)) != sizeof(float))
				bailOut("read error in file " << filename<<" vert "<<i<<","<<j);
			//if(i<20) debMsg("ReadN "<<i<<","<<j<<" "<<v[j] );
		}
		if(useNormals) normals.append(MFloatVector(v[0], v[1], v[2]));
	}
	if(!useNormals) debMsg("Not using normals!");

	// read triangles
	if (gzread(gzf, &numTris, 4) != 4)
		bailOut("read error in file " << filename);

	// store triangle indices, note - only needed for per vert UVs!
	std::vector<int> triPs;

	int addedTris = 0;
	debMsg("Reading " << numTris << " triangles");
	for (int i=0; i<numTris; i++) {
		int t[3];
		for (int j=0; j<3; j++) {
			int v;
			if (gzread(gzf, &v, 4) != 4)
				bailOut("read error in file " << filename<<" tri "<<j);
			t[j] = v;
		}
        addedTris++;
        faceCounts.append(3);
        for(int j=0; j<3; j++) {
            faceConnects.append(t[j]);
            triPs.push_back(t[j]); // store
        }        
	}
	debMsg("Added " << addedTris << " triangles");

	MObject meshObject = meshFn.create(numVerts, addedTris,
									points, faceCounts, faceConnects,
									outData, &stat);
	meshFn.setObject(meshObject);
	if(useNormals) meshFn.setNormals(normals, MSpace::kObject);

	// normal mesh loading done, now add optional parts...
	// free some memory
	points.clear();
	faceCounts.clear();
	faceConnects.clear();

	// read optional attributes
	bool fileEnd = false;
	bool haveUvs = false;
	bool haveVertcols = false;
    std::vector<float> vcR, vcG;
	
	while (!fileEnd) {
		int idIn[4] = {-1,-1,-1,-1};
		for (int j=0; j<4; j++) {
			if (gzread(gzf, &idIn[j], sizeof(int)) != 4) {
				// failed
				fileEnd= true;
			} else {
				// ok...
			} 
		}

		int type=0;
		if( idIn[0]==0 &&
		    idIn[1]=='u' &&
		    idIn[2]=='v' &&
		    idIn[3]=='w' ) type = 1; // texcoords
		else if( idIn[0]==0 &&
		    idIn[1]=='v' &&
		    idIn[2]=='c' &&
		    idIn[3]=='o' ) type = 2; // per vertex colors
		else if( idIn[0]==0 &&
            idIn[1]=='U' &&
            idIn[2]=='V' &&
            idIn[3]=='W' ) type = 3; // per vertex UVs
        else if( idIn[0]==0 &&
            idIn[1]=='v' &&
            idIn[2]=='d' &&
            idIn[3]=='e' ) type = 4; // per vertex smoke densities
        else if( idIn[0]==0 &&
            idIn[1]=='v' &&
            idIn[2]=='x' &&
            idIn[3]=='f' ) type = 5; // per vertex flags
        else if(!fileEnd) {
			//debMsg("unkown id tag "<<idIn[0]<<" "<<idIn[1]<<" "<<idIn[2]<<" "<<idIn[3]<<" " );
		}

		if (type==1) {
			debMsg("Loading tex coords...");
			haveUvs = true;

			MDagPath meshDag;
			MDagPath::getAPathTo( meshObject, meshDag );
			meshDag.extendToShape();

			// Array to hold our UVs
			MFloatArray aUVu;
			MFloatArray aUVv;

			// Array to hold UVSet names
			MStringArray        UVSetNames;
			UVSetNames.append( "map1" );
			//UVSetNames.append( "map2new" ); // test add new one...
			//meshFn.createUVSet( UVSetNames[0] ); //
			meshFn.clearUVs(); // necessary?
			int uvCnt = 0;

			for (int tri=0; tri<numTris; tri++) {
				float uvw[9];
				for (int point=0; point<3; point++) {
					for (int coord=0; coord<3; coord++) {

						if (gzread(gzf, &uvw[point*3+coord], sizeof(float)) != sizeof(float)) {
							debMsg("read error in file " << filename<<" reading texcoord "<<tri<<":"<<point<<","<<coord );
							haveUvs = false; 
						} else {
							//if(coord!=2) uvVec.push_back( uvw[point*3+coord] );
						}

					} // coord 
				}

				// set in maya
				for (int point=0; point<3; point++) { 
					// store in UV table
					meshFn.setUV(uvCnt, uvw[point*3+0], uvw[point*3+1] ); //, &UVSetNames[0] );
					uvCnt++;
				} 
				// assign UVs to polygon
				meshFn.assignUV(tri, 0, tri*3+0);
				meshFn.assignUV(tri, 1, tri*3+1);
				meshFn.assignUV(tri, 2, tri*3+2);

			} // tri
			if(haveUvs){ debMsg( "UV coords loaded, #"<<meshFn.numUVs() <<" "); }
			else { debMsg( "Error loading UVs...");}

		} else if (type==3) {
			debMsg("Loading per vert tex coords...");
			bool havePervertUVs = true;

			MDagPath meshDag;
			MDagPath::getAPathTo( meshObject, meshDag );
			meshDag.extendToShape();

			// Array to hold our UVs
			MFloatArray aUVu;
			MFloatArray aUVv;

			// Array to hold UVSet names
			MStringArray        UVSetNames;
			UVSetNames.append( "map1" );
			meshFn.clearUVs(); // necessary?
			//int uvCnt = 0;

			for (int point=0; point<numVerts; point++) {
				float uvw[2];

				for (int coord=0; coord<2; coord++) {
					if (gzread(gzf, &uvw[coord], sizeof(float)) != sizeof(float)) {
						debMsg("read error in file " << filename<<" reading texcoord "<<point<<","<<coord );
						havePervertUVs = false; 
					} 
				} // coord 

				// set in maya
				meshFn.setUV(point, uvw[0], uvw[1] ); //, &UVSetNames[0] );

			} // tri
			for (int tri=0; tri<numTris; tri++) {
				for (int point=0; point<3; point++) {
					meshFn.assignUV(tri, point, triPs[tri*3+point] );
				}
			}
			if(havePervertUVs){ debMsg( "per vert UV coords loaded, #"<<meshFn.numUVs() <<" "); }
			else { debMsg( "Error loading per vert UVs...");}

		} else if (type==2) {
            debMsg("Loading vert cols...");
            haveVertcols = true;

            MDagPath meshDag;
            MDagPath::getAPathTo( meshObject, meshDag );
            meshDag.extendToShape();
            // meshFn has to be shape!  
            MString csName("colorSet1");
            meshFn.createColorSetWithName(csName);

            for (int point=0; point<numVerts; point++) { 
                float col[3];
                for(int j=0; j<3; j++) {
                    if (gzread(gzf, &col[j], sizeof(float)) != sizeof(float)) {
                        debMsg("read error in file " << filename<<" reading vertcol "<<point<<","<<j );
                        haveVertcols = false; 
                    }  else {
                        //vcolVec.push_back( col ); 
                    }
                } // j

                MColor vcol = MColor( col[0], col[1], col[2] );
                meshFn.setVertexColor(vcol, point);
            }

            if(haveVertcols) { debMsg( "VCOL vertex color sets #"<<meshFn.numColorSets() <<" "); } 
            else {debMsg( "Error loading vertCols...");}
        } else if (type==4) {
            debMsg("Loading smoke densities in R channel...");
            haveVertcols = true;

            vcR.resize(numVerts);
            if (gzread(gzf, &vcR[0], numVerts*sizeof(float)) != sizeof(float)*numVerts) { 
                debMsg("read error in file " << filename<<" reading smoke density " );
                haveVertcols = false; 
            }
        } else if (type==5) {
            debMsg("Loading vertex flags in G channel...");
            haveVertcols = true;

            vcG.resize(numVerts);
            if (gzread(gzf, &vcG[0], numVerts*sizeof(float)) != sizeof(float)*numVerts) { 
                debMsg("read error in file " << filename<<" reading smoke density " );
                haveVertcols = false; 
            }            
        } else {
			//debMsg("No additional attributes to load...\n");
		}
	}
    if (haveVertcols) {
        MDagPath meshDag;
        MDagPath::getAPathTo( meshObject, meshDag );
        meshDag.extendToShape();
        // meshFn has to be shape!  
        MString csName("colorSet1");
        meshFn.createColorSetWithName(csName);
        meshFn.setCurrentColorSetName(csName);
        //meshFn.clearColors(&csName);
            
        for (int point=0; point<numVerts; point++) { 
            float r=0,g=0,b=0;
            if (point < vcR.size()) r=vcR[point];
            if (point < vcG.size()) g=vcG[point];
            MColor vcol = MColor(r,g,b);
            meshFn.setVertexColor(vcol, point);
        }
    }
	gzclose(gzf);


	delete filename;
	return true;
}

MStatus bobjFluidObject::compute(const MPlug& plug, MDataBlock& data)

{
	MStatus returnStatus;

	if (plug == outputMesh) {
		/* Get input data */
		MDataHandle timeData = data.inputValue( time, &returnStatus ); 
		McheckErr(returnStatus, "Error getting time data handle\n");
		MTime time = timeData.asTime();

		MDataHandle inFileMaskData = data.inputValue( inFileMask, &returnStatus ); 
		McheckErr(returnStatus, "Error getting time data handle\n");
		const char* inFileMaskString = inFileMaskData.asString().asChar();

		MDataHandle indexOffsetData = data.inputValue( indexOffset, &returnStatus ); 
		McheckErr(returnStatus, "Error getting indexOffset data handle\n");
		int indexOffsetInt = indexOffsetData.asInt();

		MDataHandle scaleData = data.inputValue( scale, &returnStatus ); 
		McheckErr(returnStatus, "Error getting scale data handle\n");
		float scaleFloat = scaleData.asFloat();

		/* Get output object */

		MDataHandle outputHandle = data.outputValue(outputMesh, &returnStatus);
		McheckErr(returnStatus, "ERROR getting polygon data handle\n");

		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);
		McheckErr(returnStatus, "ERROR creating outputData");

		// try to load current object
		if ( loadMeshData(time, inFileMaskString, indexOffsetInt, scaleFloat, newOutputData, returnStatus)) {
			McheckErr(returnStatus, "ERROR creating output mesh");
		} else {
			// loading failed, create empty mesh
			createNullMesh(newOutputData, returnStatus);
			McheckErr(returnStatus, "ERROR creating null mesh");
		}

		outputHandle.set(newOutputData);
		data.setClean( plug );
	} else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}

MStatus initializePlugin(MObject obj)
{
	MStatus   status;
	MFnPlugin plugin(obj, PYTHON_COMPANY, "3.0", "Any");

	status = plugin.registerNode("bobjFluidObject", bobjFluidObject::id,
						 bobjFluidObject::creator, bobjFluidObject::initialize);
	if (!status) {
		status.perror("registerNode");
		return status;
	}

	return status;
}

MStatus uninitializePlugin(MObject obj)
{
	MStatus	  status;
	MFnPlugin plugin(obj);

	status = plugin.deregisterNode(bobjFluidObject::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

	return status;
}
