#include "match_features.h"
#include "geometry_queries.h"

#include "read_obj.h"

#include <iostream>
#include <set>
#include <list>
#include <cfloat>

bool inverts_incident_tets(int vIndex,
                           const Vec3f& newPos,
                           const TetMesh& mesh,
                           const std::map<int,Vec3f>& movedVertices
                               = std::map<int,Vec3f>())
{
    // Check if snapping the given vertex to the given position 
    // would make any tets degenerate or inverted.
    // Include the new positions of vertices provided in movedVertices.

    std::vector<int> incidentTets;
    mesh.getIncidentTetIndices(vIndex, incidentTets);

    for (size_t i = 0; i < incidentTets.size(); ++i)
    {
        int t = incidentTets[i];
        Tet tet = mesh.getTet(t);

        // Replace any vertices that have been moved.
        for (int u = 0; u < 4; ++u)
        {
            int v = mesh.T(t)[u];
            if (v == vIndex)
            {
                tet[u] = newPos;
            }
            else
            {
                std::map<int,Vec3f>::const_iterator m = movedVertices.find(v);
                if (m != movedVertices.end())
                {
                    tet[u] == m->second;
                }
            }
        }

        // Check orientation of tet.
        int o = tet.orientation();
        if (o <= 0)
        {
            // Snapping the vertex inverts this tet!
            return true;
        }
    }
    return false;
}


// Find all boundary vertices incident to the given boundary vertex.
void
find_boundary_neighbours(int vert_idx,
                         const std::vector<int> &boundaryVerts,
                         const std::vector<Vec3i> &boundaryTris,
                         std::vector<int> &neighbours)
{
    neighbours.clear();

    std::set<int> neighbourSet;
    for (size_t t = 0; t < boundaryTris.size(); ++t) {
        if (boundaryTris[t][0] == vert_idx ||
            boundaryTris[t][1] == vert_idx ||
            boundaryTris[t][2] == vert_idx) {
            neighbourSet.insert(boundaryTris[t][0]);
            neighbourSet.insert(boundaryTris[t][1]);
            neighbourSet.insert(boundaryTris[t][2]);
        }
    }
    neighbourSet.erase(vert_idx);

    for (std::set<int>::iterator itr = neighbourSet.begin(); 
            itr != neighbourSet.end(); ++itr) {
        neighbours.push_back(*itr);
    }
}


struct SearchNode {
    int vertex;
    float pathCost;

    SearchNode(int vert, float cost) : vertex(vert), pathCost(cost) {}

    // Used to sort by cost.
    bool operator < (const SearchNode & other) const {
        return pathCost < other.pathCost;
    }
};
std::ostream& operator << (std::ostream & out, const SearchNode & node) {
    return out << "SearchNode { " << node.vertex 
        << " : pathCost = " << node.pathCost << " } ";
}


// Find an optimal path through the boundary mesh between two vertices.
// Optimality is determined by how closely the path follows the given feature.
bool
boundary_edge_search(int v1, int v2,
                     const TetMesh& mesh,
                     float dx,
                     const std::vector<int>& boundaryVerts,
                     const std::vector<Vec3i>& boundaryTris,
                     const Feature& feature,
                     const std::vector<int>& featureEndpoints,
                     const std::map<int,int>& vertexFeatureMap,
                     std::vector<int>& path,
                     std::vector<Vec3f>& projectedPositions,
                     float& pathCost,
                     bool unsafe = false)
{
    // Treat the boundary mesh like a graph whose edge weights are the 
    // distance from the vertex to which they point to the feature.
    // We use Dijkstra's algorithm to find the "shortest" path, 
    // ie: the path with the smallest edge weights.

    std::set<SearchNode> searchFront; // Nodes visited but not finished
    std::set<int> finishedNodes;      // Nodes we've finished with
    std::map<int, float> weight;      // Weights of vertices
    std::map<int, Vec3f> snappedPos;  // Snapped positions of vertices
    std::map<int, float> minPathCost; // Minimum path costs of vertices
    std::map<int, int> previous;      // Previous vertex for min-cost path
    std::vector<int> neighbours;      // Utility vector for storing neighbours.

    // Add vertices that have already been constrained to a feature to
    // our finishedNodes set.  This will prevent them from being used again.
    for (size_t i = 0; i < featureEndpoints.size(); ++i)
    {
        // Don't exclude the endpoints of the current feature. 
        if (featureEndpoints[i] != v1 && featureEndpoints[i] != v2)
        {
            finishedNodes.insert(featureEndpoints[i]);
        }
    }
    for (std::map<int,int>::const_iterator itr = vertexFeatureMap.begin();
         itr != vertexFeatureMap.end(); ++itr)
    {
        finishedNodes.insert(itr->first);
    }

    // Add the starting search node to the search front.
    SearchNode startNode(v1, 0.0f);
    assert(dist2(feature.projectToFeature(mesh.V(v1)), mesh.V(v1)) == 0.0f);
    assert(dist2(feature.projectToFeature(mesh.V(v2)), mesh.V(v2)) == 0.0f);
    weight[v1] = 0.0f;
    snappedPos[v1] = mesh.V(v1);
    snappedPos[v2] = mesh.V(v2);
    minPathCost[v1] = 0.0f;
    searchFront.insert(startNode);

    // Begin search.
    // Terminate when we find the target vertex v2, or when we have
    // fully explored v1's connected component.
    while (!searchFront.empty())
    {
        // Pop the shortest-distance node off the search front queue.
        SearchNode u = *(searchFront.begin());

        // Check if we've reached our goal.
        if (u.vertex == v2)
        {
            // We've finished searching, so recover the path between
            // v1 and v2 using the previous map.
            std::list<int> pathList;
            int prev = v2;
            while (prev != v1)
            {
                pathList.push_front(prev);
                assert(previous.find(prev) != previous.end());
                prev = previous[prev];
            }
            pathList.push_front(v1);
            
            // Construct path, projectedPositions (output parameters)
            path.clear();
            for (std::list<int>::iterator listItr = pathList.begin(); 
                 listItr != pathList.end(); ++listItr)
            {
                path.push_back(*listItr);
                projectedPositions.push_back(snappedPos[*listItr]);
            }
            pathCost = u.pathCost;
            return true;
        }

        searchFront.erase(searchFront.begin());
        finishedNodes.insert(u.vertex);

        // Visit each neighbour of u.
        find_boundary_neighbours(u.vertex, 
                                 boundaryVerts, boundaryTris, 
                                 neighbours);
        for (size_t n = 0; n < neighbours.size(); ++n)
        {
            // If the current neighbour was already visited, skip it.
            if (finishedNodes.find(neighbours[n]) != finishedNodes.end())
            {
                continue;
            }

            // Create a SearchNode to represent the neighbour, v.
            SearchNode v(neighbours[n], FLT_MAX);
            
            // Find the weight of v.
            // Stored in weight[v.vertex] if already computed.
            // Also retrieve/compute the snapped position of vertex v.
            float w;
            Vec3f p;
            std::map<int,float>::const_iterator wItr = weight.find(v.vertex);
            if (wItr == weight.end())
            {
                // Compute the weight of v. 
                // The weight is the distance squared between the position 
                // of v and the closest point to v on the feature,
                // plus a small constant factor to encourage the algorithm to 
                // choose paths with fewer edges.
                p = feature.projectToFeature(mesh.V(v.vertex));
                w = dist2(p, mesh.V(v.vertex));
                w += 1e-3; // Discourage unneccesarily long paths.
                weight[v.vertex] = w;
                snappedPos[v.vertex] = p;
            }
            else
            {
                w = wItr->second; // Weight was already computed, retrieve it.
                assert(snappedPos.find(v.vertex) != snappedPos.end());
                p = snappedPos[v.vertex];
            }

            // If moving this vertex would invert any tets,
            // disallow it from being part of the path.
            // Pass in a map of moved vertices along the path to v.
            // Only do this check if unsafe operations are not allowed.
            if (!unsafe)
            {
                std::map<int,Vec3f> movedVertices;
                int prev = u.vertex;
                while (prev != v1)
                {
                    assert(snappedPos.find(prev) != snappedPos.end());
                    assert(previous.find(prev) != previous.end());
                    movedVertices[prev] = snappedPos[prev];
                    prev = previous[prev];
                }
                if (inverts_incident_tets(v.vertex, p, mesh, movedVertices))
                {
                    continue; // Skip v and move to the next neighbour vertex.
                }
            }
            else
            {
                // Limit the total distance that a vert can move.
                if ((w-1e-3) > 16.0*dx*dx)
                {
                    continue;
                }
            }
            

            // Find the current minimum path cost to reach v.
            // Will be stored in minPathCost[v.vertex], 
            // if it has already been computed.
            std::map<int,float>::iterator costItr = minPathCost.find(v.vertex);
            if (costItr != minPathCost.end())
            {
                v.pathCost = costItr->second;
            } // Otherwise, leave it at +infinity

            // Compute the cost of reaching v from u.
            float newPathCost = u.pathCost + w;
            if (newPathCost < v.pathCost)
            {
                // Reaching v from u is a shorter path than our previous path
                // to reach v.  Update our data to reflect this.
                
                // Remove v from searchFront so we can update its pathCost.
                searchFront.erase(v);

                // Update the cost of reaching v.
                v.pathCost = newPathCost;
                if (costItr == minPathCost.end())
                {
                    minPathCost[v.vertex] = v.pathCost;
                }
                else
                {
                    costItr->second = v.pathCost;
                }

                // Update v's "previous" node to be u.
                previous[v.vertex] = u.vertex;

                // Insert the updated v back into the searchFront.
                searchFront.insert(v);
            } 
        } // for neighbours
    } // while !searchFront.empty()

    // If we got this far, then v2 is not in the connected component of v1.
    return false;
}


int
find_closest_vert(const Vec3f& point,
                  const TetMesh& mesh,
                  const std::vector<int>& boundaryVerts,
                  float dx,
                  const std::map<int,int>& alreadyUsed,
                  float& distanceSquared)
{
    int closestVert = -1;
    distanceSquared = FLT_MAX;
    float dx2 = dx*dx;
    for (size_t b = 0; b < boundaryVerts.size(); ++b)
    {
        // Find the closest boundary vertex to the given point.
        float currDist2 = dist2(point, mesh.V(boundaryVerts[b]));
        if (currDist2 < distanceSquared)
        {
            // Check that this vertex has not already been used.
            // Accept it anyway if its within dx of the point.
            if (currDist2 > dx2 && 
                alreadyUsed.find(boundaryVerts[b]) != alreadyUsed.end())
            {
                continue;
            }

            // Check that snapping this vertex will not invert a tet.
            if (inverts_incident_tets(boundaryVerts[b], point, mesh))
            {
                continue;
            }

            closestVert = boundaryVerts[b];
            distanceSquared = currDist2;
        }
    }
    return closestVert;
}


bool
match_feature_endpoints(TetMesh& mesh,
                        FeatureSet& featureSet,
                        float dx,
                        const std::vector<int>& boundaryVerts,
                        std::vector<int>& featureEndpoints,
                        std::map<Vec3f,int>& endpointVertexMap)
{
    assert(boundaryVerts.size() > 0);

    // This output variable will keep track of which boundary vertices get
    // snapped to feature endpoints.
    featureEndpoints.clear();

    // List of feature endpoints we need to match.
    std::vector<Vec3f> endpoints;
    featureSet.getFeatureEndpoints(endpoints);

    // List of endpoints sortable by distance to vertex that will be snapped.
    std::list< std::pair<float,int> > toSnap;
    std::vector<int> closestVerts; // Closest vertex to each endpoint
    closestVerts.reserve(endpoints.size());

    // Map of vertices to the endpoints they have been snapped to.
    std::map<int,int> snappedVertMap;


    //
    // Find closest vertex to each endpoint.
    //
    for (size_t ep = 0; ep < endpoints.size(); ++ep)
    {
        float distance2;
        int closestVert = find_closest_vert(endpoints[ep], mesh, boundaryVerts,
                                            dx, snappedVertMap, distance2);
        if (closestVert == -1)
        {
            // We couldn't find a vertex to snap to the endpoint that wouldn't
            // invert any tets.  This is a fatal error.
            std::cout << "ERROR: No vertex could be snapped to the feature "
                      << "endpoint at {" << endpoints[ep] << "}" << std::endl;
            return false;
        }

        closestVerts.push_back(closestVert);
        toSnap.push_back(std::make_pair(distance2, ep));
    }
    
    toSnap.sort(); // Sort endpoints by distance to vertex


    //
    // Snap vertices to endpoints in closest-first order.
    // Ensure endpoints are merged whenever there isn't a 1-to-1 mapping of
    // vertices and endpoints.
    // Also make sure we don't invert any tets when we snap.
    //
    while (!toSnap.empty())
    {
        int currEndpoint = toSnap.front().second;
        toSnap.pop_front();

        int vert = closestVerts[currEndpoint];

        // Check if the vert has already been used by a different endpoint.
        std::map<int,int>::iterator ep = snappedVertMap.find(vert);
        if (ep != snappedVertMap.end())
        {
            // The vert has already been snapped to an endpoint.
            // Merge that endpoint with the current one.
            // Merging may create new endpoints, so check for that and
            // add them to endpoints and toSnap if so.
            std::vector<Vec3f> newEndpoints;
            featureSet.mergeEndpoints(endpoints[ep->second], 
                                      endpoints[currEndpoint],
                                      newEndpoints);
            std::cout << "WARNING: Merged endpoint at {" <<
                endpoints[currEndpoint] << "} to endpoint at {" <<
                endpoints[ep->second] << "}" << std::endl;
                
            // Process newly created endpoints.
            bool newEndpointsAdded = newEndpoints.size() > 0;
            for (size_t newEp = 0; newEp < newEndpoints.size(); ++newEp)
            {
                endpoints.push_back(newEndpoints[newEp]);

                float distance2;
                int closestVert = find_closest_vert(newEndpoints[newEp], 
                                                    mesh, boundaryVerts, dx,
                                                    snappedVertMap, distance2);
                if (closestVert == -1)
                {
                    // We couldn't find a vertex to snap to the endpoint that 
                    // wouldn't invert any tets.  This is a fatal error.
                    std::cout << "ERROR: No vertex could be snapped to the " 
                              << "endpoint created in the merge: {"
                              << newEndpoints[newEp] << "}" << std::endl;
                    return false;
                }
                
                closestVerts.push_back(closestVert);
                toSnap.push_back(std::make_pair(distance2,endpoints.size()-1));
            }
            if (newEndpointsAdded)
            {
                toSnap.sort();
            }

            continue;
        }

        // Make sure snapping this vert will not invert any tets.
        if (inverts_incident_tets(vert, endpoints[currEndpoint], mesh))
        {
            // Need to find an alternative vertex that hasn't been used yet
            // and won't invert any tets.
            float distance2;
            int closestVert = find_closest_vert(endpoints[currEndpoint], mesh, 
                                                boundaryVerts, dx,
                                                snappedVertMap, distance2);
            if (closestVert == -1)
            {
                // Failed to find a valid alternative vertex! Fatal error.
                std::cout << "ERROR: No alternative vertex could be snapped "
                          << "to the feature endpoint at {" 
                          << endpoints[currEndpoint] << "}" << std::endl;
                return false;
            }

            closestVerts[currEndpoint] = closestVert;
            toSnap.push_back(std::make_pair(distance2, currEndpoint));
            toSnap.sort(); // re-sort the endpoints
            continue;
        }

        // Do the snap.
        mesh.V(vert) = endpoints[currEndpoint];
        endpointVertexMap[endpoints[currEndpoint]] = vert;
        featureEndpoints.push_back(vert);
        snappedVertMap[vert] = currEndpoint;
    }

    return true;
}


bool
match_features(TetMesh& mesh,
               FeatureSet& featureSet,
               float dx,
               const std::vector<int>& boundaryVerts,
               const std::vector<Vec3i>& boundaryTris,
               std::vector<int>& featureEndpoints,
               std::map<int, int>& vertexFeatureMap,
               bool unsafe)
{
    // Start by consolidating the feature set based on the minimum resolution
    // that we're able to resolve.
    // TODO: Make this work!
    static const double DIST_THRESHOLD_FACTOR = 0.3;
    featureSet.consolidate(DIST_THRESHOLD_FACTOR*dx);
   
    // Record a map of snapped boundary vertex locations to their indices.
    std::map<Vec3f, int> endpointVertexMap;


    //
    // First we need to match the endpoints of all the features.
    //
    bool result = match_feature_endpoints(mesh, featureSet, dx, boundaryVerts,
                                          featureEndpoints, endpointVertexMap);
    if (!result)
    {
        return false;
    }


    //
    // Find optimal paths through the boundary mesh between feature endpoints.
    //

    // The list of feature we will try to resolve.
    std::vector<Feature> features;
    featureSet.getFeatures(features);

    // Store the shortest paths we find for each feature here.
    std::vector< std::vector<int> > featurePaths(features.size());
    std::vector< std::vector<Vec3f> > snappedPositions(features.size());

    // Store the cost of each path we find.
    std::vector<float> pathCosts(features.size());

    for (size_t f = 0; f < features.size(); ++f)
    {
        const Feature &feature = features[f];

        // Features that are points should have already been handled during
        // the endpoint-snapping operation.
        if (feature.isPoint())
        {
            continue;
        }

        // Find the boundary vertices that were snapped to the current
        // feature's endpoints.
        Vec3f endpoint1 = feature.getEndpoint1();
        Vec3f endpoint2 = feature.getEndpoint2();
        
        assert(endpointVertexMap.find(endpoint1) != endpointVertexMap.end());
        assert(endpointVertexMap.find(endpoint2) != endpointVertexMap.end());
        assert(endpoint1 != endpoint2);
        assert(vertexFeatureMap.empty());
        
        int e1 = endpointVertexMap[endpoint1];
        int e2 = endpointVertexMap[endpoint2];

        // Our goal is to find a path through the boundary mesh
        // between the two feature endpoint vertices.
        // Optimally, the path should result in as little tet deformation as
        // possible when each vertex on the path is snapped to the feature.
        if (!boundary_edge_search(e1, e2, mesh, dx,
                                  boundaryVerts, boundaryTris,
                                  feature, featureEndpoints, vertexFeatureMap,
                                  featurePaths[f], snappedPositions[f],
                                  pathCosts[f]))
        {
            // Couldn't find a path between the two feature endpoints...
            // Try again with unsafe operations if allowed.
            bool success = false;
            if (unsafe)
            {
                success = boundary_edge_search(e1, e2, mesh, dx,
                                               boundaryVerts, boundaryTris,
                                               feature, featureEndpoints, 
                                               vertexFeatureMap,
                                               featurePaths[f], 
                                               snappedPositions[f],
                                               pathCosts[f],
                                               true);
            }

            if (!success)
            {
                // Ideally this shouldn't happen, but just in case:
                std::cout << "ERROR: No valid path between feature endpoints:" 
                    << std::endl << "[(" << mesh.V(e1) << "), ("  
                                     << mesh.V(e2) << ")]" << std::endl;
                featurePaths[f].clear();
            }
            else
            {
                // Path was found successfully (but might not be safe).
                assert(featurePaths[f].front() == e1);
                assert(featurePaths[f].back() == e2);
            }
        }
        else
        {
            // Path was found successfully.
            assert(featurePaths[f].front() == e1);
            assert(featurePaths[f].back() == e2);
        }
    }


    //
    // Snap vertices along each path to their corresponding feature.
    // As we do this, we need to make sure we don't snap a vertex to multiple
    // features.  Do this by snaping the "best" paths first, and re-searching
    // for any paths that use vertices that have already been snapped.
    // If we find a feature that cannot be resolved without using vertices
    // snapped to other features, only snap the unique vertices.
    // We also need to make sure that no tets are inverted when snapping.
    //
    
    // First build a sorted list of paths and their costs 
    // (paths referenced by integer index, paired with float cost).
    std::list< std::pair<float, int> > pathsToSnap;
    for (size_t p = 0; p < featurePaths.size(); ++p)
    {
        if (featurePaths[p].size() > 0) // Don't include empty paths.
        {
            pathsToSnap.push_back(std::make_pair(pathCosts[p], p));
        }
    }
    pathsToSnap.sort();

    // Keep track of paths that can still be fulfilled with unique vertices,
    // and without inverting any tets.
    std::vector<bool> unique(features.size(), true);
    std::vector<bool> preserveOrient(features.size(), true);

    // Iterate over paths until we've handled them all and removed them from
    // the pathsToSnap list.
    while (!pathsToSnap.empty())
    {
        // Extract the smallest-cost path from pathsToSnap.
        int curr = pathsToSnap.front().second;
        pathsToSnap.pop_front();

        // Snap all unique path vertices onto the current feature.
        // Skip vertices that have already been snapped to other features.
        for (size_t p = 1; p < featurePaths[curr].size()-1; ++p)
        {
            int v = featurePaths[curr][p];
            Vec3f newPos = snappedPositions[curr][p];

            if (vertexFeatureMap.find(v) != vertexFeatureMap.end())
            {
                continue; // Vertex already used, skip it.
            }

            // Make sure moving this vertex will not invert any tets.
            if (!unsafe && inverts_incident_tets(v, newPos, mesh))
            {
                continue; // This would invert tets, skip it.
            }

            // Project vertex onto the feature.
            mesh.V(v) = newPos;

            // Record which vertices were snapped to this feature.
            vertexFeatureMap[v] = curr;
        }

        // DEBUGGING
        std::vector<Vec3f> objVerts;
        std::vector<Vec3i> objTris;
        mesh.getBoundary(objVerts, objTris);
        write_objfile(objVerts, objTris, "45_feature_matching.obj");

        // Now that we've snapped vertices to features, we need to re-find
        // any paths that used the vertices that are now spoken for.
        // We also need to re-find paths that now invert tets.
        bool reSort = false;

        // Iterate over all paths left to snap, checking if they're valid.
        for (std::list< std::pair<float,int> >::iterator itr = 
             pathsToSnap.begin(); itr != pathsToSnap.end(); ++itr)
        {
            int path = itr->second;

            if (!unique[path] || !preserveOrient[path])
            {
                // We already know no valid path exists for this feature.
                continue;
            }
            
            // Check if the current path uses any vertices already in use. 
            bool uniquePath = true;
            for (size_t p = 1; p < featurePaths[path].size()-1; ++p)
            {
                int v = featurePaths[path][p];
                if (vertexFeatureMap.find(v) != vertexFeatureMap.end())
                {
                    uniquePath = false;
                    break;
                }
            }

            // Check if current path preserves orientation
            bool preservesOrientation = true;
            // Don't bother checking if we already failed uniqueness test.
            // Also don't bother if unsafe operations are allowed.
            if (uniquePath || unsafe) 
            {
                // Build map of vertices moved by this path.
                std::map<int,Vec3f> movedVertices;
                for (size_t p = 1; p < featurePaths[path].size()-1; ++p)
                {
                    int v = featurePaths[path][p];
                    movedVertices[v] = snappedPositions[path][p];
                }

                // Check each vertex for orientation preservation.
                for (size_t p = 1; p < featurePaths[path].size()-1; ++p)
                {
                    int v = featurePaths[path][p];
                    if (inverts_incident_tets(v, snappedPositions[path][p],
                                              mesh, movedVertices))
                    {
                        preservesOrientation = false;
                        break;
                    }
                }
            }

            if (!uniquePath || !preservesOrientation)
            {
                // The path uses an already-snapped vertex.
                // Try finding a new path.
                std::vector<int> altPath;
                std::vector<Vec3f> altPositions;
                float altCost;
               
                assert(featurePaths[path].front() == 
                    endpointVertexMap[features[path].getEndpoint1()]);
                assert(featurePaths[path].back() == 
                    endpointVertexMap[features[path].getEndpoint2()]);

                // featurePaths[path].front()/.back() are the endpoints of
                // features[path] expressed as boundary mesh indices,
                // which is what boundary_edge_search expects as input.
                if (boundary_edge_search(featurePaths[path].front(),
                                         featurePaths[path].back(),
                                         mesh, dx,
                                         boundaryVerts, boundaryTris,
                                         features[path], 
                                         featureEndpoints, vertexFeatureMap,
                                         altPath, altPositions, altCost))
                {
                    // We found an alternate path, save it and update the
                    // current entry in pathsToSnap (pointed to by itr).
                    featurePaths[path] = altPath;
                    snappedPositions[path] = altPositions;
                    pathCosts[path] = altCost;
                    itr->first = altCost;
                    reSort = true;
                }
                else
                {
                    bool success = false;
                    if (unsafe)
                    {
                        success = boundary_edge_search(
                                                    featurePaths[path].front(),
                                                    featurePaths[path].back(),
                                                    mesh, dx,
                                                    boundaryVerts, boundaryTris,
                                                    features[path], 
                                                    featureEndpoints, 
                                                    vertexFeatureMap,
                                                    altPath, altPositions, 
                                                    altCost,
                                                    true);
                        if (success)
                        {
                            featurePaths[path] = altPath;
                            snappedPositions[path] = altPositions;
                            pathCosts[path] = altCost;
                            itr->first = altCost;
                            reSort = true;
                        }
                    }

                    if (!success)
                    {
                        // Could not find a valid alternate path.
                        unique[path] = uniquePath;
                        preserveOrient[path] = preservesOrientation;

                        std::cout << "WARNING: Unable to fully match feature ";
                        if (!uniquePath)
                        {
                            std::cout << "(vertices already used) ";
                        }
                        else
                        {
                            std::cout << "(tets inverted) ";
                        }
                        std::cout << "with endpoints:" << std::endl
                            << "[(" << mesh.V(featurePaths[path].front()) 
                            << "), (" << mesh.V(featurePaths[path].back()) 
                            << ")]" << std::endl;
                    }
                }
            }

        } // end for pathsToSnap.begin()

        // We may have modified the costs of the contents of pathsToSnap,
        // so we need to re-sort it.
        if (reSort)
        {
            pathsToSnap.sort();
        }

    } // end while !pathsToSnap.empty()

    return true;
}


