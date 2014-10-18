#ifndef TET_MESH_H
#define TET_MESH_H

#include "vec.h"
#include <vector>

//
// Data structure representing a single tetrahedron.
//
struct Tet
{
    Vec3f v[4]; // The four vertices of the tetrahedron.

    Tet() { }
    Tet(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const Vec3f& v3);

    Vec3f& operator[] (int i) { return v[i]; }
    const Vec3f& operator[] (int i) const { return v[i]; }

    int orientation(double &sixSignedVolume) const;
    int orientation() const;

    float volume() const;
    float aspectRatio() const;

    void dihedralAngles(std::vector<float>& angles) const;
    void faceAngles(std::vector<float>& angles) const;

    float maxEdgeLength() const;
};


//
// Data structure representing a tetrahedral mesh.
//
struct TetMesh
{
    // Default constructor
    TetMesh() : dirty(false) {}
    // Copy constructor
    TetMesh(const TetMesh& mesh) : v(mesh.v), t(mesh.t), dirty(true) {}
    // Construct from vertices and tet indices
    TetMesh(const std::vector<Vec3f>& verts, 
            const std::vector<Vec4i>& tets) :
        v(verts), t(tets), dirty(true) {}

    // Returns the position of the vertex stored at index i.
    Vec3f& V(int i) { return v[i]; }
    const Vec3f& V(int i) const { return v[i]; }

    // Returns the indices of the tetrahedron stored at index i.
    Vec4i& T(int i) { dirty = true; return t[i]; }
    const Vec4i& T(int i) const { return t[i]; }

    // Returns the list of this tet mesh's vertices.
    std::vector<Vec3f>& verts() { dirty = true; return v; }
    const std::vector<Vec3f>& verts() const { return v; }

    // Returns the list of this tet mesh's tetrahedron indices.
    std::vector<Vec4i>& tets() { dirty = true; return t; }
    const std::vector<Vec4i>& tets() const { return t; }

    // Returns the number of vertices in the tet mesh.
    size_t vSize() const { return v.size(); }

    // Returns the number of tetrahedra in the mesh.
    size_t tSize() const { return t.size(); }

    // Returns the tetrahedron stored at index t.
    Tet getTet(int t) const;

    // Fills the given list with all indices of all the vertices 
    // that are adjacent to the given vertex.
    void getAdjacentVertices(int vertex, std::vector<int>& vertices) const;
    // Fills the given list with all the tets incident to the given vertex.
    void getIncidentTets(int vertex, std::vector<Tet>& tets) const;
    // Fills the given list with the indices of tets incident to vertex.
    void getIncidentTetIndices(int vertex, std::vector<int>& tetIndices) const;

    // Clean up the mesh, getting rid of unused vertices.
    // Tetrahedral indices adjusted accordingly.
    void compactMesh();

    // Extract boundary vertices and triangles from the tet mesh.
    // Vertices will be returned in boundary_verts as indices into the 
    // tet mesh vertex list v.
    // Triangles are stored in boundary_tris.  These are also indices into v, 
    // NOT as indices into boundary_verts!!!
    // Assumes that the tet mesh is well-formed and closed...
    void getBoundary(std::vector<int>& boundary_verts,
                     std::vector<Vec3i>& boundary_tris) const;

    // Extract boundary vertices and triangles from the tet mesh.
    // Vertices will be copied into boundary_verts.
    // Triangles are stored in boundary_tris.  These are indices into 
    // boundary_verts. 
    // Assumes that the tet mesh is well-formed and closed...
    void getBoundary(std::vector<Vec3f>& boundary_verts,
                     std::vector<Vec3i>& boundary_tris) const;

    // Write this TetMesh to file with the given name.
    // The output file will be of the .TET format.
    bool writeToFile(const char* filename) const;

    // Compute information about this TetMesh (dihedral angles, tet volumes,
    // etc.) and write it to the file with the given name.
    bool writeInfoToFile(const char* filename) const;

private:
    std::vector<Vec3f> v;
    std::vector<Vec4i> t;

    // Stores the indices of the tetrahedra incident to each vertex.
    mutable std::vector< std::vector<int> > incidenceMap;
    // Indicates whether the incidence map may need to be rebuilt.
    mutable bool dirty; 

    // Re-create the incidence map from the current v and t lists.
    // This is const because we need to be able to call it from 
    // getIncidentTets(), which is const.  
    // incidenceMap and dirty are mutable as a result.
    void buildIncidenceMap() const;
};

#endif
