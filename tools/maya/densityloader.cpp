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

#include "zlib.h"
#include <maya/MIOStream.h>
#include <maya/MObject.h>
#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MFnPlugin.h>
#include <maya/MDataBlock.h>
#include <maya/MFnFluid.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MTime.h>
#include <maya/MSelectionList.h>

class FluidGridObject : public MPxNode
{
public:
					FluidGridObject() {};
	virtual 		~FluidGridObject() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject inTime, outTime;
	static MObject inFileMask;
	static MObject scale;
	static MObject indexOffset;
	static MObject fluidName;
	static MTypeId id;

protected:
	void loadDensity(const char* mask, int offset, MTime& time, MObject& obj);
};


MObject FluidGridObject::inTime;
MObject FluidGridObject::outTime;
MObject FluidGridObject::inFileMask;
MObject FluidGridObject::indexOffset;
MObject FluidGridObject::fluidName;
MObject FluidGridObject::scale;
MTypeId FluidGridObject::id( 0x80012 );

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg << endl;				\
		return MS::kFailure;		\
	}

void* FluidGridObject::creator()
{
	return new FluidGridObject;
}

MStatus FluidGridObject::initialize()
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
	FluidGridObject::inTime = unitAttr.create( "inTime", "itm", MFnUnitAttribute::kTime, 0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating FluidGridObject time attribute\n");
	returnStatus = addAttribute(FluidGridObject::inTime);
	McheckErr(returnStatus, "ERROR adding time attribute\n");

	FluidGridObject::outTime = unitAttr.create( "outTime", "otm", MFnUnitAttribute::kTime, 0.0, &returnStatus );
	McheckErr(returnStatus, "ERROR creating FluidGridObject time attribute\n");
	unitAttr.setWritable(false);
	unitAttr.setStorable(false);
	returnStatus = addAttribute(FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR adding time attribute\n");
	
	string = stringData.create(&returnStatus);
	McheckErr(returnStatus, "ERROR creating string\n");
	FluidGridObject::inFileMask = typedAttr.create( "inFileMask", "fil", MFnData::kString, string, &returnStatus);
	McheckErr(returnStatus, "ERROR creating FluidGridObject inFileMask attribute\n");
	returnStatus = addAttribute(FluidGridObject::inFileMask);
	McheckErr(returnStatus, "ERROR adding inFileMask attribute\n");

	FluidGridObject::indexOffset = numericAttr.create( "indexOffset", "ofs", MFnNumericData::kInt, -1, &returnStatus);
	McheckErr(returnStatus, "ERROR creating FluidGridObject indexOffset attribute\n");
	returnStatus = addAttribute(FluidGridObject::indexOffset);
	McheckErr(returnStatus, "ERROR adding indexOffset attribute\n");

	FluidGridObject::scale = numericAttr.create( "scale", "sc", MFnNumericData::kFloat, 1.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating FluidGridObject scale attribute\n");
	returnStatus = addAttribute(FluidGridObject::scale);
	McheckErr(returnStatus, "ERROR adding scale attribute\n");

	string = stringData.create(&returnStatus);
	McheckErr(returnStatus, "ERROR creating string\n");
	FluidGridObject::fluidName = typedAttr.create( "fluidName", "fluid", MFnData::kString, string, &returnStatus);	
	McheckErr(returnStatus, "ERROR creating FluidGridObject fluidname attribute\n");
	returnStatus = addAttribute(FluidGridObject::fluidName);
	McheckErr(returnStatus, "ERROR adding fluid attribute\n");
		
	/*
	 * Set dependencies
	 */	
	returnStatus = attributeAffects(FluidGridObject::inTime, FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR in making attribute time affect output\n");
	returnStatus = attributeAffects(FluidGridObject::scale, FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR in making attribute scale affect output\n");
	returnStatus = attributeAffects(FluidGridObject::inFileMask, FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR in making attribute filemask affect output\n");
	returnStatus = attributeAffects(FluidGridObject::indexOffset, FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR in making attribute indexOffset affect output\n");
	returnStatus = attributeAffects(FluidGridObject::fluidName, FluidGridObject::outTime);
	McheckErr(returnStatus, "ERROR in making attribute fluidname affect output\n");
	
	return MS::kSuccess;
}

MStatus FluidGridObject::compute(const MPlug& plug, MDataBlock& data) {
	MStatus returnStatus;
	if (plug == outTime) {
		/* Get input data */
		MDataHandle timeData = data.inputValue( inTime, &returnStatus ); 
		McheckErr(returnStatus, "Error getting time data handle\n");
		MTime time = timeData.asTime();

		MDataHandle inFileMaskData = data.inputValue( inFileMask, &returnStatus ); 
		McheckErr(returnStatus, "Error getting filemask handle\n");
		const char* inFileMaskString = inFileMaskData.asString().asChar();

		MDataHandle fluidNameData = data.inputValue( fluidName, &returnStatus ); 
		McheckErr(returnStatus, "Error getting fluidname handle\n");
		const char* fluidNameString = fluidNameData.asString().asChar();

		MDataHandle indexOffsetData = data.inputValue( indexOffset, &returnStatus ); 
		McheckErr(returnStatus, "Error getting indexOffset data handle\n");
		int indexOffsetInt = indexOffsetData.asInt();

		MDataHandle scaleData = data.inputValue( scale, &returnStatus ); 
		McheckErr(returnStatus, "Error getting scale data handle\n");
		float scaleFloat = scaleData.asFloat();

		/* Get output object */
		MDataHandle outputHandle = data.outputValue(outTime, &returnStatus);
		McheckErr(returnStatus, "ERROR getting outTime data handle\n");
		outputHandle.set(time);
		data.setClean(plug);	

		// Obtain fluid shape object from name
		MSelectionList tempList;
		tempList.add( fluidNameString );
		if ( tempList.length() > 0 ) {
			MObject obj;
			tempList.getDependNode( 0, obj );
			loadDensity(inFileMaskString, indexOffsetInt, time, obj);
		}
			
	} else
		return MS::kUnknownParameter;
	
	return MS::kSuccess;
} 

typedef struct {
	int dimX, dimY, dimZ;
	int frames, elements, elementType, bytesPerElement, bytesPerFrame;
} DHeader;

typedef struct {
    int dimX, dimY, dimZ;
    int gridType, elementType, bytesPerElement;
} UniHeader;


void FluidGridObject::loadDensity(const char* mask, int offset, MTime& time, MObject& obj) {
	int frameNo;
	char* filename = NULL;

	// compute frame index
	frameNo = (int)time.as( MTime::kFilm );
	frameNo += offset;

	// validate input file mask
	if (!mask || strlen(mask)<3) return;
	
	const char* perc = strrchr(mask, '%');
	if (perc == NULL) return; // invalid mask string
	if (!((perc[1]=='d')||(perc[2]=='d')||(perc[3]=='d'))) return; // invalid mask string

	// get input file name
	size_t namelen = strlen(mask)+40;
	filename = new char[namelen];
#ifdef _WIN32
	sprintf_s(filename, namelen, mask, frameNo);
#else
    snprintf(filename, namelen, mask, frameNo);
#endif

	// obtain reference to Fluid Object
	if (obj.apiType() != MFn::kFluid) {
        cerr << "can't obtain fluid node handle"<< endl;
	    return; 
    }
        
	MFnFluid fluid(obj);
		
	// read file
	gzFile gzf = gzopen(filename, "rb");
	if (!gzf) {
        cerr << "can't open file "<< filename << endl;	    
        return;
    }
	cerr << "Loading file " << filename << endl;
    
	// Header
    int sx, sy, sz;
    float* field=0;
	char id[5]; 
	gzread(gzf, &id, 4); id[4]='\0';
	if (!strcmp(id, "DDF2")) {        
        // old DDF format
        DHeader hdr;
        gzread(gzf, &hdr, sizeof(DHeader));
        if (hdr.elements != 1 || hdr.elementType != 1 || hdr.bytesPerElement != 4) { // only 1-el float fields supported
            gzclose(gzf);
            return; 
        }
        // load flags
        sx = hdr.dimX; sy = hdr.dimY; sz = hdr.dimZ;
        int fieldsize = sx*sy*sz;
        int sizeFlag = fieldsize * sizeof(char);
        char* flags = (char*)malloc(sizeFlag);
        int numRead = gzread(gzf, flags, sizeFlag);
        if (numRead != sizeFlag) {
            gzclose(gzf);
            free(flags);
            return;
        }
        
        // load values
        int sizeField = fieldsize*sizeof(float);
        field = (float*)malloc(sizeField);
        numRead = gzread(gzf, field, sizeField);
        gzclose(gzf);
        if (numRead != sizeField) {
            free(flags);
            free(field);
            return;
        }
        cerr << "DDF2 format, read ok" << endl;    
    } else if(!strcmp(id, "MNT1")) {
        UniHeader hdr;
        gzread(gzf, &hdr, sizeof(UniHeader));
        if (hdr.elementType != 1) { gzclose(gzf); return; }
        
        sx = hdr.dimX; sy = hdr.dimY; sz = hdr.dimZ;
        int sizeField = sx*sy*sz*sizeof(float);
        field = (float*)malloc(sizeField);
        int numRead = gzread(gzf, field, sizeField);
        gzclose(gzf);
        if (numRead != sizeField) {
            free(field); return;
        }
        cerr << "MNT1 format, read ok" << endl;    
    }

	// set dimensions
	unsigned xres, yres, zres;
	fluid.getResolution(xres, yres, zres);
	if (sx != xres || sy != yres || sz != zres) {
		const double asz = 5.0;
		int mindim = (sx < sy) ? sx : sy;
		if (sz < mindim) mindim = sz;
		double mult = asz / (double) mindim;
		fluid.setSize(sx, sy, sz, mult*(double)sx, mult*(double)sy, mult*(double)sz, false);
	}
	// set properties
	fluid.setVelocityMode(MFnFluid::kZero, MFnFluid::kConstant);
	fluid.setFuelMode(MFnFluid::kZero, MFnFluid::kConstant);
	fluid.setTemperatureMode(MFnFluid::kZero, MFnFluid::kConstant);
	fluid.setDensityMode(MFnFluid::kStaticGrid, MFnFluid::kConstant);
	
	// copy values
	float* dest = fluid.density();
	if (dest != NULL) {
		unsigned totalSize =fluid.gridSize();
		if (sx*sy*sz == totalSize) {
			for (int k=0, idx=0; k<sz; k++)
				for (int j=0; j<sy; j++)
					for (int i=0; i<sx; i++) {
						dest[fluid.index(i,j,k)] = field[idx++];
					}
			fluid.updateGrid();
		}		
	}
	free(field);			
}

MStatus initializePlugin(MObject obj)
{
	MStatus   status;
	MFnPlugin plugin(obj, PYTHON_COMPANY, "6.0", "Any");

	status = plugin.registerNode("fluidGridObject", FluidGridObject::id, FluidGridObject::creator, FluidGridObject::initialize);
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

	status = plugin.deregisterNode(FluidGridObject::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}
	
	return status;
}
