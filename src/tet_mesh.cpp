#include "tet_mesh.h"
#include "geometry_queries.h"
#include <set>
#include <map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cfloat>

//
// Tet
//

Tet::Tet(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2, const Vec3f& v3)
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
}


int Tet::orientation(double &sixSignedVolume) const
{
    return orientation3D(v[0][0], v[0][1], v[0][2],
                         v[1][0], v[1][1], v[1][2],
                         v[2][0], v[2][1], v[2][2],
                         v[3][0], v[3][1], v[3][2],
                         sixSignedVolume);
}


int Tet::orientation() const
{
    double sixSignedVolume;
    return orientation(sixSignedVolume);
}


float Tet::volume() const
{
    double sixSignedVolume;
    orientation(sixSignedVolume);
    return sixSignedVolume/6.0;
}


float Tet::aspectRatio() const
{
    // Find the two edges with the largest cross product magnitude.
    // ie: max(cross(w,x))
    Vec3f edges[6];
    edges[0] = v[1] - v[0];
    edges[1] = v[2] - v[0];
    edges[2] = v[3] - v[0];
    edges[3] = v[2] - v[1];
    edges[4] = v[3] - v[1];
    edges[5] = v[3] - v[2];

    float maxCrossMag2 = 0.0;
    for (int i = 0; i < 5; ++i)
    {
        for (int j = i+1; j < 6; ++j)
        {
            float crossMag2 = mag2(cross(edges[i], edges[j]));
            if (maxCrossMag2 < crossMag2)
            {
                maxCrossMag2 = crossMag2;
            }
        }
    }
    float maxCrossMag = sqrt(maxCrossMag2);
    
    // We also need six times the volume of the tet, which we can
    // get from the orientation primitive.
    double sixSignedVolume;
    orientation(sixSignedVolume);

    // Last thing we need is the maximum edge length.
    float maxEdge = maxEdgeLength();

    // aspectRatio = (l_max * max(cross(w,x)) / 6*V
    return (maxEdge * maxCrossMag) / sixSignedVolume;
}


void Tet::dihedralAngles(std::vector<float>& angles) const
{
    angles.clear();
    angles.resize(6);

    // Indices for the vertices of each face of the tet.
    static const Vec3i FACE_INDICES[4] = {
        Vec3i(0, 1, 2),
        Vec3i(0, 2, 3),
        Vec3i(0, 3, 1), 
        Vec3i(1, 3, 2)
    };
    bool degenerate[4];

    // First compute the normal of each face.
    Vec3f normals[4];
    for (int f = 0; f < 4; ++f)
    {
        Vec3f e1 = v[FACE_INDICES[f][1]] - v[FACE_INDICES[f][0]];
        Vec3f e2 = v[FACE_INDICES[f][2]] - v[FACE_INDICES[f][0]];
        normals[f] = cross(e1, e2);
        float m = mag(normals[f]);
        if (m > 0.0)
        {
            normals[f] /= m;
            degenerate[f] = false;
        }
        else
        {
            degenerate[f] = true;
        }
    }

    // Then compute the angle between each pair of normals.
    int a = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i+1; j < 4; ++j)
        {
            if (degenerate[i] || degenerate[j])
            {
                angles[a] = 0.0;
            }
            else
            {
                float cosAngle = -dot(normals[i], normals[j]);
                angles[a] = acos(cosAngle);
            }
            ++a;
        }
    }
}


void Tet::faceAngles(std::vector<float>& angles) const
{
    angles.clear();
    angles.resize(12);

    // For each vertex, compute the three face angles at that vertex.
    int a = 0;
    for (int vert = 0; vert < 4; ++vert)
    {
        Vec3f e1 = v[(vert+1)%4] - v[vert];
        Vec3f e2 = v[(vert+2)%4] - v[vert];
        Vec3f e3 = v[(vert+3)%4] - v[vert];
        float m1 = mag(e1);
        float m2 = mag(e2);
        float m3 = mag(e3);
        if (m1 > 0.0) e1 /= m1;
        if (m2 > 0.0) e2 /= m2;
        if (m3 > 0.0) e3 /= m3;

        if (m1 > 0.0 && m2 > 0.0)
            angles[a++] = acos(dot(e1, e2));
        else
            angles[a++] = 0.0;
        if (m2 > 0.0 && m3 > 0.0)
            angles[a++] = acos(dot(e2, e3));
        else
            angles[a++] = 0.0;
        if (m3 > 0.0 && m1 > 0.0)
            angles[a++] = acos(dot(e3, e1));
        else
            angles[a++] = 0.0;
    }
}


float Tet::maxEdgeLength() const
{
    float maxLength2 = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = i+1; j < 4; ++j)
        {
            float edgeLength2 = dist2(v[i], v[j]);
            if (maxLength2 < edgeLength2)
            {
                maxLength2 = edgeLength2;
            }
        }
    }
    return sqrt(maxLength2);
}



//
// TetMesh
//

Tet TetMesh::getTet(int i) const
{
    const Vec4i& tet = t[i];
    return Tet(v[tet[0]],
               v[tet[1]],
               v[tet[2]],
               v[tet[3]]);
}


void TetMesh::getAdjacentVertices(int vertex, std::vector<int>& vertices) const
{
    if (dirty || incidenceMap.size() == 0)
    {
        buildIncidenceMap();
    }

    assert(incidenceMap.size() == v.size());
    assert(vertex < static_cast<int>(incidenceMap.size()));

    vertices.clear();

    const std::vector<int>& incidentTets = incidenceMap[vertex];
    for (size_t i = 0; i < incidentTets.size(); ++i)
    {
        for (int u = 0; u < 4; ++u)
        {
            // Add the current vertex if it hasn't been added already.
            int vIndex = t[incidentTets[i]][u];
            bool alreadyAdded = false;
            for (size_t j = 0; j < vertices.size(); ++j)
            {
                if (vertices[j] == vIndex)
                {
                    alreadyAdded = true;
                    break;
                }
            }
            if (!alreadyAdded && vIndex != vertex)
            {
                vertices.push_back(vIndex);
            }
        }
    }
}


void TetMesh::getIncidentTets(int vertex, std::vector<Tet>& tets) const
{
    if (dirty || incidenceMap.size() == 0)
    {
        buildIncidenceMap();
    }

    assert(incidenceMap.size() == v.size());
    assert(vertex < static_cast<int>(incidenceMap.size()));

    tets.clear();

    const std::vector<int>& incidentTets = incidenceMap[vertex];
    for (size_t i = 0; i < incidentTets.size(); ++i)
    {
        tets.push_back(getTet(incidentTets[i]));
    }
}


void TetMesh::getIncidentTetIndices(int vertex, 
                                    std::vector<int>& tetIndices) const
{
    if (dirty || incidenceMap.size() == 0)
    {
        buildIncidenceMap();
    }

    assert(incidenceMap.size() == v.size());
    assert(vertex < static_cast<int>(incidenceMap.size()));

    tetIndices = incidenceMap[vertex];
}


void TetMesh::compactMesh()
{
    // Go through our tetrahedra, relabeling vertices.
    std::vector<int> remap(v.size(), -1);
    int nv = 0;
    for (size_t tet = 0; tet < t.size(); ++tet)
    {
        for (int u = 0; u < 4; ++u)
        {
            int i = t[tet][u];
            if (remap[i] == -1)
            {
                remap[i] = nv;
                ++nv;
            }
            t[tet][u] = remap[i];
        }
    }

    // Now go through the vertices, remapping them using their new labels.
    // Any vertices not found in tetrahedra will not be included.
    std::vector<Vec3f> vnew(nv);
    for (size_t i = 0; i < v.size(); ++i)
    {
        if (remap[i] >= 0) vnew[remap[i]] = v[i];
    }
    v.swap(vnew);

    dirty = true;
}


void TetMesh::getBoundary(std::vector<int>& boundary_verts,
                          std::vector<Vec3i>& boundary_tris) const
{
    boundary_verts.clear();
    boundary_tris.clear();

    // Go through every triangle of every tetrahedron.
    // If the same triangle has already been found
    // it is not a boundary face.
    std::set<Vec3i> boundary_set;
    Vec3i tet_tris[4];
    Vec3i permuted_tris[3];

    for (size_t tet = 0; tet < t.size(); ++tet)
    {
        // Consider each triangle
        tet_tris[0] = Vec3i(t[tet][0], t[tet][1], t[tet][2]);
        tet_tris[1] = Vec3i(t[tet][0], t[tet][2], t[tet][3]);
		tet_tris[2] = Vec3i(t[tet][0], t[tet][3], t[tet][1]);
        tet_tris[3] = Vec3i(t[tet][1], t[tet][3], t[tet][2]);
	
        // If the winding is wrong on the boundary, try this instead
        //tet_tris[0] = Vec3i t1(t[tet][0], t[tet][2], t[tet][1]);
		//tet_tris[1] = Vec3i t2(t[tet][0], t[tet][3], t[tet][2]);
		//tet_tris[2] = Vec3i t3(t[tet][0], t[tet][1], t[tet][3]);
        //tet_tris[3] = Vec3i t4(t[tet][1], t[tet][2], t[tet][3]);
		
		for (int tri = 0; tri < 4; ++tri)
        {
            // If the current triangle already exists in our boundary_set,
            // it will have opposite winding and an arbitrary first vertex.
            // Check all possible valid permutations.
            permuted_tris[0] = Vec3i(tet_tris[tri][0], 
                                     tet_tris[tri][2],
                                     tet_tris[tri][1]);
            permuted_tris[1] = Vec3i(tet_tris[tri][1], 
                                     tet_tris[tri][0],
                                     tet_tris[tri][2]);
            permuted_tris[2] = Vec3i(tet_tris[tri][2], 
                                     tet_tris[tri][1],
                                     tet_tris[tri][0]);
            // Attempt to erase each permutation from the set.
            // If erasure was not successful, it means the triangle isn't
            // in the set yet, so add it.
            if (!boundary_set.erase(permuted_tris[0]) &&
                !boundary_set.erase(permuted_tris[1]) &&
                !boundary_set.erase(permuted_tris[2]))
            {
                // Add the triangle to the set.
                boundary_set.insert(tet_tris[tri]);
            }
        }
    }

    // Once we're done traversing the mesh, we should have all boundary
    // triangles in boundary_set.
    // Build the set of boundary vertices, and store the boundary triangles
    // and boundary vertices in boundary_verts and boundary_tris.
    std::set<int> vertex_set;
    for (std::set<Vec3i>::iterator itr = boundary_set.begin();
         itr != boundary_set.end(); ++itr)
    {
        boundary_tris.push_back(*itr);
        
        vertex_set.insert((*itr)[0]);
        vertex_set.insert((*itr)[1]);
        vertex_set.insert((*itr)[2]);
    }
    for (std::set<int>::iterator itr = vertex_set.begin();
            itr != vertex_set.end(); ++itr)
    {
        boundary_verts.push_back(*itr);
    }
}


void TetMesh::getBoundary(std::vector<Vec3f>& boundary_verts,
                          std::vector<Vec3i>& boundary_tris) const
{
    boundary_verts.clear();
    boundary_tris.clear();

    std::vector<int> vertIndices;
    std::vector<Vec3i> triIndices;
    getBoundary(vertIndices, triIndices);

    boundary_verts.reserve(vertIndices.size());
    boundary_tris.reserve(triIndices.size());

    std::map<int, int> meshToBoundaryMap;

    for (size_t vIndex = 0; vIndex < vertIndices.size(); ++vIndex)
    {
        boundary_verts.push_back(v[vertIndices[vIndex]]);
        meshToBoundaryMap[vertIndices[vIndex]] = vIndex;
    }
    for (size_t tIndex = 0; tIndex < triIndices.size(); ++tIndex)
    {
        Vec3i tri;
        for (int i = 0; i < 3; ++i)
        {
            assert(meshToBoundaryMap.find(triIndices[tIndex][i]) !=
                   meshToBoundaryMap.end());
            tri[i] = meshToBoundaryMap[triIndices[tIndex][i]];
        }
        boundary_tris.push_back(tri);
    }
}


bool TetMesh::writeToFile(const char* filename) const
{
    FILE *out=std::fopen(filename, "w");
    if (!out) return false;
    fprintf(out, "tet %d %d\n", (int)v.size(), (int)t.size());
    for (size_t i=0; i<v.size(); ++i)
        fprintf(out, "%.7g %.7g %.7g\n", 
                v[i][0], v[i][1], v[i][2]);
    for (size_t tet=0; tet<t.size(); ++tet)
        fprintf(out, "%d %d %d %d\n", 
                t[tet][0], t[tet][1], t[tet][2], t[tet][3]);
    std::fclose(out);
    return true;
}


bool TetMesh::writeInfoToFile(const char* filename) const
{
    std::ofstream out(filename);
    if (!out.good()) return false;

    //
    // Compute dihedral angle data.
    //
    static const int DIHEDRAL_NUM_BINS = 18;
    static const float DIHEDRAL_BIN_MAX[DIHEDRAL_NUM_BINS] = 
        { 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 110.0,
          120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 175.0, 180.0};
    
    // Accumulate dihedral angles
    std::vector<float> dihedralAngles;
    dihedralAngles.reserve(t.size()*6);
    for (size_t tIdx = 0; tIdx < t.size(); ++tIdx)
    {
        Tet currTet = getTet(tIdx);
        std::vector<float> currAngles;
        currTet.dihedralAngles(currAngles);
        for (size_t i = 0; i < currAngles.size(); ++i)
        {
            dihedralAngles.push_back(currAngles[i] * 180.0 / M_PI);
        }
    }
    
    // Process dihedral angles, generating histogram and min/max.
    float minDihedral = dihedralAngles[0];
    float maxDihedral = minDihedral;
    std::vector<int> dihedralHistogram(DIHEDRAL_NUM_BINS, 0);
    for (size_t i = 0; i < dihedralAngles.size(); ++i)
    {
        float currAngle = dihedralAngles[i];
        if (currAngle < minDihedral) minDihedral = currAngle;
        if (currAngle > maxDihedral) maxDihedral = currAngle;
        for (int bin = 0; bin < DIHEDRAL_NUM_BINS; ++bin)
        {
            if (currAngle <= DIHEDRAL_BIN_MAX[bin])
            {
                ++(dihedralHistogram[bin]);
                break;
            }
        }
    }


    //
    // Compute volume data.
    //

    // Compute minimum and maximum volume.
    float minVolume = FLT_MAX;
    float maxVolume = -FLT_MAX;
    for (size_t tIdx = 0; tIdx < t.size(); ++tIdx)
    {
        Tet currTet = getTet(tIdx);
        float currVolume = currTet.volume();
        if (currVolume < minVolume) minVolume = currVolume;
        if (currVolume > maxVolume) maxVolume = currVolume;
    }


    //
    // Compute aspect ratio data.
    //
    static const int ASPECT_NUM_BINS = 12;
    static const float ASPECT_BIN_MAX[ASPECT_NUM_BINS] = 
        { 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 15.0,
          25.0, 50.0, 100.0, FLT_MAX};
    
    // Accumulate aspect ratio
    std::vector<float> aspectRatios;
    aspectRatios.reserve(t.size());
    for (size_t tIdx = 0; tIdx < t.size(); ++tIdx)
    {
        Tet currTet = getTet(tIdx);
        aspectRatios.push_back(currTet.aspectRatio());    
    }

    // Process aspect ratios, generating histogram and min/max.
    float minAspect = aspectRatios[0];
    float maxAspect = minAspect;
    std::vector<int> aspectHistogram(ASPECT_NUM_BINS, 0);
    for (size_t i = 0; i < aspectRatios.size(); ++i)
    {
        float currAspect = aspectRatios[i];
        if (currAspect < minAspect) minAspect = currAspect;
        if (currAspect > maxAspect) maxAspect = currAspect;
        for (int bin = 0; bin < ASPECT_NUM_BINS; ++bin)
        {
            if (currAspect <= ASPECT_BIN_MAX[bin])
            {
                ++(aspectHistogram[bin]);
                break;
            }
        }
    }

    
    //
    // Write out data
    //

    out << "Tetrahedral Mesh Quality Information" << std::endl << std::endl;

    out << "Mesh vertices: " << v.size() << std::endl;
    out << "Mesh tetrahedra: " << t.size() << std::endl;
    out << std::endl;

    out << std::setw(24) << std::left << "Smallest volume:" << 
           std::setprecision(5) << std::setw(8) << std::right << minVolume 
        << "   |   "
        << std::setw(24) << std::left << "Largest volume:" <<
           std::setprecision(5) << std::setw(8) << std::right << maxVolume 
        << std::endl;
    out << std::setw(24) << std::left << "Smallest aspect ratio:" << 
           std::setprecision(5) << std::setw(8) << std::right << minAspect 
        << "   |   "
        << std::setw(24) << std::left << "Largest aspect ratio:" <<
           std::setprecision(5) << std::setw(8) << std::right << maxAspect 
        << std::endl;
    out << std::setw(24) << std::left << "Smallest dihedral:" << 
           std::setprecision(5) << std::setw(8) << std::right << minDihedral 
        << "   |   "
        << std::setw(24) << std::left << "Largest dihedral:" <<
           std::setprecision(5) << std::setw(8) << std::right << maxDihedral 
        << std::endl;
    out << std::endl;

    out << "Aspect ratio histogram:" << std::endl;
    size_t halfAspectIdx = aspectHistogram.size() / 2 + 
                           aspectHistogram.size() % 2;
    for (size_t i = 0; i < halfAspectIdx; ++i)
    {
        float min = 0.0;
        if (i > 0)
        {
            min = ASPECT_BIN_MAX[i-1];
        }
        float max = ASPECT_BIN_MAX[i];

        out << std::setw(6) << std::setprecision(5) << std::right << min
            << " - "
            << std::setw(11) << std::setprecision(5) << std::left << max
            << ":" << std::setw(11) << std::right << aspectHistogram[i];
        
        out << "   |   ";
        
        size_t rightIdx = i+halfAspectIdx;
        if (rightIdx >= aspectHistogram.size())
        {
            out << std::endl;
            continue;
        }
        min = ASPECT_BIN_MAX[rightIdx-1];
        max = ASPECT_BIN_MAX[rightIdx];
        
        out << std::setw(6) << std::setprecision(5) << std::right << min
            << " - "
            << std::setw(11) << std::setprecision(5) << std::left << max
            << ":" << std::setw(11) << std::right << aspectHistogram[rightIdx];
        out << std::endl;
    }
    out << "(Aspect ratio is longest edge divided by shortest altitude)" 
        << std::endl << std::endl;

    out << "Dihedral angle histogram:" << std::endl;
    size_t halfDihedralIdx = dihedralHistogram.size() / 2 + 
                             dihedralHistogram.size() % 2;
    for (size_t i = 0; i < halfDihedralIdx; ++i)
    {
        float min = 0.0;
        if (i > 0)
        {
            min = DIHEDRAL_BIN_MAX[i-1];
        }
        float max = DIHEDRAL_BIN_MAX[i];

        out << std::setw(6) << std::setprecision(3) << std::right << min
            << " - "
            << std::setw(3) << std::right << max
            << " degrees:"
            << std::setw(11) << std::right << dihedralHistogram[i];
        
        out << "   |   ";

        size_t rightIdx = i+halfDihedralIdx;
        if (rightIdx >= dihedralHistogram.size())
        {
            out << std::endl;
            continue;
        }
        min = DIHEDRAL_BIN_MAX[rightIdx-1];
        max = DIHEDRAL_BIN_MAX[rightIdx];

        out << std::setw(6) << std::setprecision(3) << std::right << min
            << " - "
            << std::setw(3) << std::right << max
            << " degrees:"
            << std::setw(11) << std::right << dihedralHistogram[rightIdx];
        out << std::endl;
    }
    out << std::endl;

    return true;
}


void TetMesh::buildIncidenceMap() const
{
    incidenceMap.clear();
    incidenceMap.resize(v.size());

    // For each vertex, find and store all the tetrahedra that contain it.
    for (size_t vIndex = 0; vIndex < v.size(); ++vIndex)
    {
        for (size_t tIndex = 0; tIndex < t.size(); ++tIndex)
        {
            for (int i = 0; i < 4; ++i)
            {
                if (t[tIndex][i] == static_cast<int>(vIndex))
                {
                    incidenceMap[vIndex].push_back(tIndex);
                }
            }
        }
    }

    dirty = false;
}

