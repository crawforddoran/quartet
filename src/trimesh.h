#ifndef TRIMESH_H
#define TRIMESH_H

#include "vec.h"
#include <vector>


class TriMesh
{
public:
    struct HalfEdge;

    struct Vertex
    {
        Vec3f p;
        HalfEdge* e;
    };

    struct Face
    {
        HalfEdge* e;
    };

    struct HalfEdge
    {
        Vertex* v;
        Face* f;
        HalfEdge* prev, *next;
        HalfEdge* sym;
    };

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<HalfEdge> edges;

    TriMesh(const std::vector<Vec3f>& x,
            const std::vector<Vec3i>& tri);
};

#endif

