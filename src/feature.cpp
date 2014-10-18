#include "feature.h"
#include <fstream>
#include <set>
#include <list>
#include <cfloat>
#include <cstdarg>

//
// Feature
//

Vec3f Feature::projectToFeature(const Vec3f &p) const
{
    // Feature defined as sequence of vertices connected by line segments.
    // Find the closest line segment to the point, and then project onto it.
    if (vertices.size() == 0)
    {
        return p; // Nothing to do here.
    }
    else if (vertices.size() == 1)
    {
        // The feature is just a point.
        return vertices[0];
    }
    else
    {
        // The feature is comprised of one or more edges.
        // Find the nearest one and project onto that.
        Vec3f projectionPoint;
        float minDist2 = FLT_MAX;
        for (size_t i = 0; i < vertices.size()-1; ++i)
        {
            const Vec3f& e1 = vertices[i];
            const Vec3f& e2 = vertices[i+1];

            // Compute the closest point on the current edge.
            Vec3f edge(e2 - e1);
            double edgeLength2 = mag2(edge);
            float t = static_cast<float>(dot(e2-p, edge)/edgeLength2);
            // Clamp t-value between 0 and 1.
            if (t < 0.0) { t = 0.0; } else if (t > 1.0) { t = 1.0; }
            // Compute distance to closest point and save if it's smallest.
            Vec3f closestPoint = t*e1 + (1-t)*e2;
            float d2 = dist2(p, closestPoint);
            if (d2 < minDist2)
            {
                projectionPoint = closestPoint;
                minDist2 = d2;
            }
        }

        // projectionPoint now contains the closest point on any of the
        // edges to the given point.
        return projectionPoint;
    }
}


Vec3f Feature::getTangent(const Vec3f& p) const
{
    if (vertices.size() < 2)
    {
        return Vec3f(0,0,0); // No tangent exists.
    }
    else
    {
        // Find the nearest edge and return it.
        Vec3f tangentEdge(0,0,0);
        float minDist2 = FLT_MAX;
        for (size_t i = 0; i < vertices.size()-1; ++i)
        {
            const Vec3f& e1 = vertices[i];
            const Vec3f& e2 = vertices[i+1];

            Vec3f edge(e2 - e1);
            double edgeLength2 = mag2(edge);
            float t = static_cast<float>(dot(e2-p, edge)/edgeLength2);
            // Clamp t-value between 0 and 1.
            if (t < 0.0) { t = 0.0; } else if (t > 1.0) { t = 1.0; }
            // Compute distance to closest point
            Vec3f closestPoint = t*e1 + (1-t)*e2;
            float d2 = dist2(p, closestPoint);
            if (d2 < minDist2)
            {
                tangentEdge = edge;
                minDist2 = d2;
            }
        }
        assert(mag2(tangentEdge) > 0.0);

        // tangentEdge is now the edge vector of the closest edge to the point
        return normalized(tangentEdge);
    }
}


//
// FeatureSet
//

FeatureSet::FeatureSet(const std::vector<Feature>& f)
{
    addFeatures(f);
}


void FeatureSet::addFeatures(const std::vector<Feature>& features)
{
    for (size_t i = 0; i < features.size(); ++i)
    {
        addFeature(features[i]);
    }
}


void FeatureSet::addFeature(const Feature& feature)
{
    std::vector<int> featureIndices;

    for (size_t i = 0; i < feature.vertices.size(); ++i)
    {
        const Vec3f& vert = feature.vertices[i];
        int vertIndex = -1;

        // Find out if vert has already been added to vertices.
        for (size_t v = 0; v < vertices.size(); ++v)
        {
            if (vertices[v] == vert)
            {
                vertIndex = v;
            }
        }
        
        if (vertIndex == -1)
        {
            // Vert was not already added, add it.
            featureIndices.push_back(vertices.size());
            vertexMap[vert] = vertices.size();
            vertices.push_back(vert);
        }
        else
        {
            // Sanity check: make sure same vertex not added twice in a row.
            if (i > 0 && vertIndex == featureIndices[featureIndices.size()-1])
            {
                continue;
            }

            // Vert was already added, use its index.
            featureIndices.push_back(vertIndex);

            // We only want endpoints to be shared between features,
            // so split the feature being added at the current vertex.
            if (i > 0 && i < feature.vertices.size()-1)
            {
                indices.push_back(featureIndices);
                checkForLoop(indices.size()-1);
                featureIndices.clear();
                featureIndices.push_back(vertIndex);
            }

            // Also split any other features that used the current vertex.
            splitFeatures(vertIndex);
        }
    }

    if (featureIndices.size() > 0)
    {
        indices.push_back(featureIndices);
        checkForLoop(indices.size()-1);
    }
}


bool FeatureSet::readFromFile(const char *filename_format, ...)
{
    // Max line size when reading .feat files.
    static const int LINESIZE = 1024;

    // Construct filename from arguments and initialize file stream.
	va_list ap;
    va_start(ap, filename_format);
#ifdef WIN32
#define FILENAMELENGTH 256
   char filename[FILENAMELENGTH];
   _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
   std::ifstream input(filename, std::ifstream::binary);
#else
   char *filename;
   if (vasprintf(&filename, filename_format, ap) == -1) {
        return false;
   }
   std::ifstream input(filename, std::ifstream::binary);
   std::free(filename);
#endif	
    va_end(ap);

    if (!input.good()) return false;

    std::vector<Feature> features;

    // The feature file begins with the word "feat".
    // This is followed by the total number of features,
    // followed by a list of features.
    // Each individual feature begins with a single integer indicating the
    // number of vertices in that feature.
    // This is followed by a list of vertices in the feature, one per line.
    // Empty lines and lines beginning with '#' are ignored.
    char line[LINESIZE];

    // First we expect the word "feat"
    if (input.get() != 'f') return false;
    if (input.get() != 'e') return false;
    if (input.get() != 'a') return false;
    if (input.get() != 't') return false;

    // The first number is the total number of polylines.
    int numFeatures = 0;
    while (input.good()) {
        input.getline(line, LINESIZE);
        if (line[0] != 0 && line[0] != '\r' && line[0] != '#') {
            std::sscanf(line, "%d", &numFeatures);
            break;
        }
    }

    Feature feat;
    Vec3f v;
    int numVertices = 0;
    for (int f = 0; f < numFeatures && input.good(); ++f) {

        // First find the number of vertices.
        while (input.good()) {
            input.getline(line, LINESIZE);
            if (line[0] != 0 && line[0] != '\r' && line[0] != '#') {
                std::sscanf(line, "%d", &numVertices);
                break;
            }
        }

        // Then read in the vertices.
        int i = 0;
        while (i < numVertices && input.good()) {
            input.getline(line, LINESIZE);
            if (line[0] != 0 && line[0] != '\r' && line[0] != '#') {
                std::sscanf(line, "%f %f %f", &v[0], &v[1], &v[2]);
                feat.vertices.push_back(v);
                ++i;
            }
        }

        // Finally, save the feature.
        if (feat.vertices.size() > 0) {
            features.push_back(feat);
            feat.vertices.clear();
        }
    }

    // Now add the features to our set.
    addFeatures(features);

    return true;
}


bool FeatureSet::writeToFile(const char *filename_format, ...)
{
    // Construct filename from arguments and initialize file stream.
	va_list ap;
    va_start(ap, filename_format);
#ifdef WIN32
#define FILENAMELENGTH 256
   char filename[FILENAMELENGTH];
   _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
   std::ofstream output(filename, std::ifstream::binary);
#else
   char *filename;
   if (vasprintf(&filename, filename_format, ap) == -1) {
        return false;
   }
   std::ofstream output(filename, std::ifstream::binary);
   std::free(filename);
#endif	
    va_end(ap);

    if (!output.good()) return false;

    // First line is "feat" and the number of features
    output << "feat " << indices.size() << std::endl;

    // Subsequent lines are the features themselves, lists of vertices.
    for (size_t f = 0; f < indices.size(); ++f)
    {
        output << indices[f].size() << std::endl;
        for (size_t v = 0; v < indices[f].size(); ++v)
        {
            output << vertices[indices[f][v]][0] << " " 
                   << vertices[indices[f][v]][1] << " "
                   << vertices[indices[f][v]][2] << std::endl;
        }
    }

    return true;
}


void FeatureSet::consolidate(float dist)
{
    // TODO
}


void FeatureSet::getFeature(unsigned int i, Feature& feature) const
{
    feature.vertices.clear();
    feature.vertices.reserve(indices[i].size());
    for (size_t v = 0; v < indices[i].size(); ++v)
    {
        feature.vertices.push_back(vertices[indices[i][v]]);
    }
}


void FeatureSet::getFeatures(std::vector<Feature>& features) const
{
    features.clear();
    features.reserve(indices.size());

    Feature feat;
    for (size_t f = 0; f < indices.size(); ++f)
    {
        feat.vertices.clear();
        feat.vertices.reserve(indices[f].size());
        for (size_t v = 0; v < indices[f].size(); ++v)
        {
            feat.vertices.push_back(vertices[indices[f][v]]);
        }
        features.push_back(feat);
    }
}


void FeatureSet::getFeatureEndpoints(std::vector<Vec3f>& endpoints) const
{
    endpoints.clear();

    std::set<int> endpointIndices;
    int e1, e2;

    for (size_t f = 0; f < indices.size(); ++f)
    {
        assert(indices[f].size() > 0);

        // Add the first endpoint if it hasn't been already.
        e1 = indices[f][0];
        if (endpointIndices.insert(e1).second)
        {
            endpoints.push_back(vertices[e1]);
        }
        
        if (indices[f].size() > 1)
        {
            // Add the second if it hasn't been already.
            e2 = indices[f][(indices[f].size()-1)];
            if (endpointIndices.insert(e2).second)
            {
                endpoints.push_back(vertices[e2]);
            }
        }
    }
}


void FeatureSet::mergeEndpoints(const Vec3f& endpoint1, 
                                const Vec3f& endpoint2,
                                std::vector<Vec3f>& newEndpoints)
{
    newEndpoints.clear();

    assert(vertexMap.find(endpoint1) != vertexMap.end());
    assert(vertexMap.find(endpoint2) != vertexMap.end());

    int index1 = vertexMap[endpoint1];
    int index2 = vertexMap[endpoint2];

    assert(endpoint1 == vertices[index1]);
    assert(endpoint2 == vertices[index2]);

    //
    // Given two valid indices of feature endpoint vertices,
    // merge the second one to the first one.
    //

    // Replace all instances of the second one with the first one.
    // Because both are endpoints, we only need to check the endpoints of
    // features.
    for (size_t f = 0; f < indices.size(); ++f)
    {
        assert(indices[f].size() > 0);

        bool featureModified = false;
        if (indices[f][0] == index2)
        {
            indices[f][0] = index1;
            featureModified = true;
        }

        if (indices[f].size() > 1)
        {
            int endIdx = indices[f].size()-1;
            if (indices[f][endIdx] == index2)
            {
                indices[f][endIdx] = index1;
                featureModified = true;
            }
        }

        if (featureModified)
        {
            // A loop or degenerate feature may have formed, handle it.
            bool featureSplit = checkForLoop(f, dist2(endpoint1, endpoint2));
            if (featureSplit)
            {
                int endIdx = indices[f].size()-1;
                newEndpoints.push_back(vertices[indices[f][endIdx]]);
            }
        }
    }
}


void FeatureSet::autoDetectFeatures(const TriMesh& mesh, 
                                    float angle_threshold_degrees)
{
    // Go through all edges of the tri-mesh, and check if their dihedral angle
    // is smaller than angle_threshold.

    // Store half-edges whose symmetric pair has already been visited, 
    // so we don't get duplicates.
    std::set<const TriMesh::HalfEdge*> visited;
    std::vector<const TriMesh::HalfEdge*> featureEdges;

    for (size_t e = 0; e < mesh.edges.size(); ++e)
    {
        const TriMesh::HalfEdge& edge = mesh.edges[e];

        if (visited.erase(&edge))
        {
            // If this edge's symmetric pair was already visited, skip it.
            // We can erase it because we won't encounter this edge again.
            continue;
        }

        // 
        // Compute dihedral angle.
        //

        // Compute first normal
        Vec3f v1 = edge.v->p;
        Vec3f v2 = edge.next->v->p;
        Vec3f v3 = edge.prev->v->p;
        Vec3f n1 = cross(v2-v1,v3-v1);
        normalize(n1);

        // Compute second normal
        assert(edge.sym);
        v1 = edge.sym->v->p;
        v2 = edge.sym->next->v->p;
        v3 = edge.sym->prev->v->p;
        Vec3f n2 = cross(v2-v1,v3-v1);
        normalize(n2);

        // dot(n1,n2) = cos(180-theta) where theta is the dihedral angle
        // this is equivalent to cos(theta) = -dot(n1,n2)
        float dihedral_angle = acos(-dot(n1,n2)) * 180.0 / M_PI;

        if (dihedral_angle <= angle_threshold_degrees)
        {
            featureEdges.push_back(&edge);
        }

        // Store symmetric halfedge so we don't do it again later.
        visited.insert(edge.sym);
    }

    // We now have a list of all edges with angle smaller than our threshold.
    // Figure out which ones should be combined into a single feature based
    // on connectivity and angle between them.
    // Then add the features to the given feature set.
    std::vector<bool> used(featureEdges.size(), false);
    for (size_t e = 0; e < featureEdges.size(); ++e)
    {
        if (used[e]) continue; // Don't re-use edges

        std::list<const TriMesh::Vertex*> feature;
        feature.push_back(featureEdges[e]->v);
        feature.push_back(featureEdges[e]->next->v);
        used[e] = true;

        // Go through all other remaining edges and see if they are connected.
        for (size_t e2 = e+1; e2 < featureEdges.size(); ++e2)
        {
            if (used[e2]) continue; // Don't re-use edges

            const TriMesh::HalfEdge* newEdge = featureEdges[e2];

            static const float POLYLINE_ANGLE_THRESHOLD = M_PI/4.0;
            static const float POLYLINE_THRESHOLD = 
                                cos(POLYLINE_ANGLE_THRESHOLD);

            // Check connectivity to back of current feature.
            if (newEdge->v == feature.back())
            {
                // Compute angle between the newEdge and the back edge.
                Vec3f newE = newEdge->next->v->p - newEdge->v->p;
                normalize(newE);
                Vec3f backE = feature.back()->p - (*(++feature.rbegin()))->p;
                normalize(backE);
                if (dot(newE, backE) > POLYLINE_THRESHOLD)
                {
                    feature.push_back(newEdge->next->v);
                    used[e2] = true;
                    e2 = e; // The feature has changed, restart the search.
                }
            }
            else if (newEdge->next->v == feature.back())
            {
                // Compute angle between the newEdge and the back edge.
                Vec3f newE = newEdge->v->p - newEdge->next->v->p;
                normalize(newE);
                Vec3f backE = feature.back()->p - (*(++feature.rbegin()))->p;
                normalize(backE);
                if (dot(newE, backE) > POLYLINE_THRESHOLD)
                {
                    feature.push_back(newEdge->v);
                    used[e2] = true;
                    e2 = e; // The feature has changed, restart the search.
                }
            }

            // Check connectivity to front of current feature.
            else if (newEdge->v == feature.front())
            {
                // Compute angle between the newEdge and the front edge.
                Vec3f newE = newEdge->next->v->p - newEdge->v->p;
                normalize(newE);
                Vec3f frontE = feature.front()->p - (*(++feature.begin()))->p;
                normalize(frontE);
                if (dot(newE, frontE) > POLYLINE_THRESHOLD)
                {
                    feature.push_front(newEdge->next->v);
                    used[e2] = true;
                    e2 = e; // The feature has changed, restart the search.
                }
            }
            else if(newEdge->next->v == feature.front())
            {
                // Compute angle between the newEdge and the front edge.
                Vec3f newE = newEdge->v->p - newEdge->next->v->p;
                normalize(newE);
                Vec3f frontE = feature.front()->p - (*(++feature.begin()))->p;
                normalize(frontE);
                if (dot(newE, frontE) > POLYLINE_THRESHOLD)
                {
                    feature.push_front(newEdge->v);
                    used[e2] = true;
                    e2 = e; // The feature has changed, restart the search.
                }
            }
        } // end for e2 < featureEdges.size()

        // We should now have a complete feature, so add it to the feature set.
        Feature f;
        for (std::list<const TriMesh::Vertex*>::iterator itr = feature.begin();
             itr != feature.end(); ++itr)
        {
            f.vertices.push_back((*itr)->p);
        }
        addFeature(f);
    } // end for e < featureEdges.size()

}


bool FeatureSet::checkForLoop(int fIndex, float collapseRadius2)
{
    // Check if the given feature is a loop, and if so handle it.

    if (indices[fIndex].size() <= 1) {
        return false; // Feature is empty or a point, so not a loop.
    }

    int endIdx = indices[fIndex].size()-1;
    if (indices[fIndex][0] == indices[fIndex][endIdx])
    {
        // The feature is a loop.

        // Check if the feature is degenerate.
        if (indices[fIndex].size() == 2)
        {
            // Convert the feature into a point.
            indices[fIndex].pop_back();
        }
        else
        {
            const Vec3f& endpoint = vertices[indices[fIndex][0]];

            // Find the vertex on the feature furthest away from the endpoints.
            size_t splitIndex = 1;
            float splitDist2 = dist2(vertices[indices[fIndex][1]], endpoint);
            for (int v = 2; v < endIdx; ++v)
            {
                float d2 = dist2(vertices[indices[fIndex][v]], endpoint);
                if (d2 > splitDist2)
                {
                    splitIndex = v;
                    splitDist2 = d2;
                }
            }
            
            // If the furthest vertex is within collapseRadius, 
            // collapse the whole feature to a point.
            if (splitDist2 <= collapseRadius2)
            {
                while (indices[fIndex].size() > 1)
                {
                    indices[fIndex].pop_back();
                }
                return false; // No new features or endpoints added.
            }

            // Split the feature at splitIndex.
            std::vector<int> splitFeature;
            for (size_t i = splitIndex; i < indices[fIndex].size(); ++i)
            {
                splitFeature.push_back(indices[fIndex][i]);
            }
            while (indices[fIndex].size() > splitIndex+1)
            {
                indices[fIndex].pop_back();
            }
            indices.push_back(splitFeature);

            return true; // Feature was split and new feature added.
        }
    }

    return false;
}


void FeatureSet::splitFeatures(int vIndex)
{
    for (size_t f = 0; f < indices.size(); ++f)
    {
        assert(indices[f].size() > 0);
        for (size_t v = 1; v < indices[f].size()-1; ++v)
        {
            if (indices[f][v] == vIndex)
            {
                // Found the given vertex index in the interior of a feature.
                // Split the feature at v.
                std::vector<int> splitFeature;

                for (size_t i = v; i < indices[f].size(); ++i)
                {
                    splitFeature.push_back(indices[f][i]);
                }
                while (indices[f].size() > v+1)
                {
                    indices[f].pop_back();
                }
                indices.push_back(splitFeature);
            }
        }
    }
}


