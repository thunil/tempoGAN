/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Meshes
 * 
 *  note: this is only a temporary solution, details are bound to change
 *        long term goal is integration with Split&Merge code by Wojtan et al.
 *
 ******************************************************************************/

#ifndef _MESH_H
#define _MESH_H

#include <vector>
#include "manta.h"
#include "vectorbase.h"
#include <set>
namespace Manta {

// fwd decl
class GridBase;
class LevelsetGrid;
class FlagGrid;
class MACGrid;
class Shape;

//! Node position and flags
struct Node {
    Node() : flags(0), pos(Vec3::Zero), normal(Vec3::Zero) {}
    Node(const Vec3& p) : flags(0), pos(p) {}
    int flags;
    Vec3 pos, normal;
};

//! Carries indices of its nodes
struct Triangle {    
    Triangle() : flags(0) { c[0] = c[1] = c[2] = 0; }
    Triangle(int n0, int n1, int n2) : flags(0) { c[0]=n0; c[1]=n1; c[2]=n2; }
    
    int c[3];
    int flags;
};

//! For fast access to nodes and neighboring triangles
struct Corner {
    Corner() : tri(-1), node(-1), opposite(-1), next(-1), prev(-1) {};
    Corner(int t, int n) : tri(t), node(n), opposite(-1), next(-1), prev(-1) {}    
    
    int tri;
    int node;
    int opposite;
    int next;
    int prev;
};

//! Base class for mesh data channels (texture coords, vorticity, ...)
struct NodeChannel {
    virtual ~NodeChannel() {};
    virtual void resize(int num) = 0;
    virtual int size() = 0;
    virtual NodeChannel* clone() = 0;
    
    virtual void addInterpol(int a, int b, Real alpha) = 0;
    virtual void mergeWith(int node, int delnode, Real alpha) = 0;
    virtual void renumber(const std::vector<int>& newIndex, int newsize) = 0;
};

//! Node channel using only a vector
template<class T>
struct SimpleNodeChannel : public NodeChannel {
    SimpleNodeChannel() {};
    SimpleNodeChannel(const SimpleNodeChannel<T>& a) : data(a.data) {}
    void resize(int num) { data.resize(num); }
    virtual int size() { return data.size(); }
    virtual void renumber(const std::vector<int>& newIndex, int newsize);
    
    //virtual void addSplit(int from, Real alpha) { data.push_back(data[from]); }    

    std::vector<T> data;
};

//! Base class for mesh data channels (texture coords, vorticity, ...)
struct TriChannel {
    virtual ~TriChannel() {};
    virtual void resize(int num) = 0;
    virtual TriChannel* clone() = 0;
    virtual int size() = 0;
    
    virtual void addNew() = 0;
    virtual void addSplit(int from, Real alpha) = 0;
    virtual void remove(int tri) = 0;
};

//! Tri channel using only a vector
template<class T>
struct SimpleTriChannel : public TriChannel {
    SimpleTriChannel() {};
    SimpleTriChannel(const SimpleTriChannel<T>& a) : data(a.data) {}
    void resize(int num) { data.resize(num); }
    void remove(int tri) { if (tri!=(int)data.size()-1) data[tri] = *data.rbegin(); data.pop_back(); }
    virtual int size() { return data.size(); }
        
    virtual void addSplit(int from, Real alpha) { data.push_back(data[from]); }    
    virtual void addNew() { data.push_back(T()); }
    
    std::vector<T> data;
};

struct OneRing {
    OneRing() {}
    std::set<int> nodes;
    std::set<int> tris;
};

//! Triangle mesh class
/*! note: this is only a temporary solution, details are bound to change
          long term goal is integration with Split&Merge code by Wojtan et al.*/
PYTHON() class Mesh : public PbClass {
public:
    PYTHON() Mesh(FluidSolver* parent);
    virtual ~Mesh();
    virtual Mesh* clone();
    
    enum NodeFlags { NfNone = 0, NfFixed = 1, NfMarked = 2, NfKillme = 4, NfCollide = 8 };
    enum FaceFlags { FfNone = 0, FfDoubled = 1, FfMarked = 2 };
    enum MeshType { TypeNormal = 0, TypeVortexSheet };
    
    virtual MeshType getType() { return TypeNormal; }
        
    Real computeCenterOfMass(Vec3& cm) const;
    void computeVertexNormals();
    
    // plugins
    PYTHON() void clear();
    PYTHON() void load (std::string name, bool append = false);
    PYTHON() void fromShape (Shape& shape, bool append = false);
    PYTHON() void save (std::string name);
    PYTHON() void advectInGrid(FlagGrid& flags, MACGrid& vel, int integrationMode);
    PYTHON() void scale(Vec3 s);
    PYTHON() void offset(Vec3 o);

	PYTHON() void computeLevelset(LevelsetGrid& levelset, Real sigma, Real cutoff=-1.);
	//! map mesh to grid with sdf
	PYTHON() void applyMeshToGrid(GridBase* grid, FlagGrid* respectFlags=0, Real cutoff=-1.);
    
    // ops
    Mesh& operator=(const Mesh& o);
    
    // accessors
    inline int numTris() const { return mTris.size(); }
    inline int numNodes() const { return mNodes.size(); }
    inline int numTriChannels() const { return mTriChannels.size(); }
    inline int numNodeChannels() const { return mNodeChannels.size(); }
    
    inline Triangle& tris(int i) { return mTris[i]; }
    inline Node& nodes(int i) { return mNodes[i]; }    
    inline Corner& corners(int tri, int c) { return mCorners[tri*3+c]; }
    inline Corner& corners(int c) { return mCorners[c]; }
    inline NodeChannel* nodeChannel(int i) { return mNodeChannels[i]; }
    inline TriChannel* triChannel(int i) { return mTriChannels[i]; }

	// allocate memory (eg upon load)
	void resizeTris(int numTris);
	void resizeNodes(int numNodes);
    
    inline bool isNodeFixed(int n) { return mNodes[n].flags & NfFixed; }
    inline bool isTriangleFixed(int t) { return (mNodes[mTris[t].c[0]].flags & NfFixed) || (mNodes[mTris[t].c[1]].flags & NfFixed) || (mNodes[mTris[t].c[2]].flags & NfFixed); }
    
    inline const Vec3 getNode(int tri, int c) const { return mNodes[mTris[tri].c[c]].pos; }
    inline Vec3& getNode(int tri, int c) { return mNodes[mTris[tri].c[c]].pos; }
    inline const Vec3 getEdge(int tri, int e) const { return getNode(tri,(e+1)%3) - getNode(tri,e); }
    inline OneRing& get1Ring(int node) { return m1RingLookup[node]; }
    inline Real getFaceArea(int t) const { Vec3 c0 = mNodes[mTris[t].c[0]].pos; return 0.5*norm(cross(mNodes[mTris[t].c[1]].pos - c0, mNodes[mTris[t].c[2]].pos - c0)); }
    inline Vec3 getFaceNormal(int t) { Vec3 c0 = mNodes[mTris[t].c[0]].pos; return getNormalized(cross(mNodes[mTris[t].c[1]].pos - c0, mNodes[mTris[t].c[2]].pos - c0)); }
    inline Vec3 getFaceCenter(int t) const { return (mNodes[mTris[t].c[0]].pos + mNodes[mTris[t].c[1]].pos + mNodes[mTris[t].c[2]].pos) / 3.0; }
    inline std::vector<Node>& getNodeData() { return mNodes; }
    
    void mergeNode(int node, int delnode);
    int addNode(Node a);
    int addTri(Triangle a);
    void addCorner(Corner a);
    void removeTri(int tri);
    void removeTriFromLookup(int tri);
    void removeNodes(const std::vector<int>& deletedNodes);
    void rebuildCorners(int from=0, int to=-1);
    void rebuildLookup(int from=0, int to=-1);
    void rebuildQuickCheck();
    void fastNodeLookupRebuild(int corner);
    void sanityCheck(bool strict=true, std::vector<int>* deletedNodes=0, std::map<int,bool>* taintedTris=0);
    
    void addTriChannel(TriChannel* c) { mTriChannels.push_back(c); rebuildChannels(); }
    void addNodeChannel(NodeChannel* c) { mNodeChannels.push_back(c); rebuildChannels(); }

protected:    
    void rebuildChannels();
    
    std::vector<Node> mNodes;
    std::vector<Triangle> mTris;
    std::vector<Corner> mCorners;
    std::vector<NodeChannel*> mNodeChannels;
    std::vector<TriChannel*> mTriChannels;
    std::vector<OneRing> m1RingLookup;
};




// ***************************************************************************************************************
// Implementation

template<class T>
void SimpleNodeChannel<T>::renumber(const std::vector<int>& newIndex, int newsize) {
    for(size_t i=0; i<newIndex.size(); i++) {
        if(newIndex[i]!=-1)
            data[newIndex[i]] = data[newsize+i];
    }
    data.resize(newsize);
}



} //namespace
#endif
