#include "trimesh.h"

#include <map>


TriMesh::TriMesh(const std::vector<Vec3f>& x,
                 const std::vector<Vec3i>& tri)
    : vertices(x.size()), faces(tri.size()), edges(tri.size()*3)
{
    for (size_t v = 0; v < vertices.size(); ++v)
    {
        vertices[v].p = x[v];
    }

    // For keeping track of edges that we haven't paired off yet.
    std::map<std::pair<int,int>,int> unmatchedEdges;
    std::map<std::pair<int,int>,int>::iterator edgeItr;

    for (size_t f = 0; f < faces.size(); ++f)
    {
        const Vec3i& t = tri[f];
        HalfEdge* e = &(edges[3*f]);

        for (int i = 0; i < 3; ++i)
        {
            // Associate vertex and edge.
            e[i].v = &(vertices[t[i]]);
            e[i].v->e = &(e[i]);

            // Associate edge with face.
            e[i].f = &(faces[f]);

            // Link up edges.
            e[i].prev = &(e[(i+2)%3]);
            e[i].next = &(e[(i+1)%3]);

            // Match symmetric edges
            std::pair<int,int> tempEdge(t[i], t[(i+1)%3]);
            if (tempEdge.first > tempEdge.second)
            {
                // Lower vertex index always comes first.
                tempEdge.first = tempEdge.second;
                tempEdge.second = t[i];
            }
            edgeItr = unmatchedEdges.find(tempEdge);
            if (edgeItr == unmatchedEdges.end())
            {
                // The pair to this edge has not been found yet.
                // Save the index of the current half-edge for later.
                unmatchedEdges[tempEdge] = 3*f+i;
            }
            else
            {
                // The pair was already found, pair it with the current edge.
                e[i].sym = &(edges[edgeItr->second]);
                e[i].sym->sym = &(e[i]);
                unmatchedEdges.erase(edgeItr);
            }
        }

        // Associate face with the first edge.
        e[0].f->e = &(e[0]);
    }
}


