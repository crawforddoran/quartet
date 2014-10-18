#include "optimize_tet_mesh.h"

#include <iostream>

#include <set>
#include <queue>
#include <cfloat>

#include "tet_quality.h"


// Return a random number in the range -D to D
float
randf(float D)
{
    return D*(2*static_cast<float>(rand())/static_cast<float>(RAND_MAX) - 1);
}


// Do local optimization of the position of the given vertex.
float
optimize_vertex(int vertex,
                TetMesh& mesh,
				const SDF &sdf,
                const FeatureSet& featureSet,
				bool isBoundary, 
                int featureIndex)
{
	static const int NUM_RAND_PERTURBATIONS = 16;
	std::vector<Vec3f> testx;

    const float P=0.5*sdf.dx; // max perturbation in each coordinate

    if (featureIndex != -1) {
        // Optimize along the feature to which the vertex is constrained.
        // Generate random perturbation points along feature.
        // Do this by distributing points along the feature's tangent vector,
        // then projecting them all onto the feature.
        Feature feature;
        featureSet.getFeature(featureIndex, feature);
        Vec3f tangent = feature.getTangent(mesh.V(vertex));
        normalize(tangent);
        for (int i = 0; i < NUM_RAND_PERTURBATIONS; ++i) {
            testx.push_back(feature.projectToFeature(mesh.V(vertex) +
                                                     randf(P)*tangent));
        }
        // and include original point and extreme points along feature.
        testx.push_back(feature.projectToFeature(mesh.V(vertex)));
        testx.push_back(feature.projectToFeature(mesh.V(vertex)-P*tangent));
        testx.push_back(feature.projectToFeature(mesh.V(vertex)+P*tangent));

    } else if (isBoundary) {
        // Generate random perturbation points on tangent plane, in
        // a square of side length dx, then projected onto isosurface
        Vec3f normal = sdf.normal(mesh.V(vertex));
        if (mag2(normal) == 0.0)
            normal = Vec3f(0.0, 1.0, 0.0); // Default to y-axis if no normal.
        Vec3f u(1,0,0);
        if(std::fabs(normal[0])>.5)
            u = Vec3f(-normal[1], normal[0], 0.f);
        else
            u = Vec3f(0.f, -normal[2], normal[1]);
        u /= mag(u);
        Vec3f v = cross(normal, u);
        for (int i = 0; i < NUM_RAND_PERTURBATIONS; ++i) {
            Vec3f t = mesh.V(vertex) + randf(P)*u + randf(P)*v;
            testx.push_back(sdf.projectToIsosurface(t));
        }
        // also original point and projected corners of unit square
        testx.push_back(sdf.projectToIsosurface(mesh.V(vertex)));
        testx.push_back(sdf.projectToIsosurface(mesh.V(vertex)-P*u-P*v));
        testx.push_back(sdf.projectToIsosurface(mesh.V(vertex)+P*u-P*v));
        testx.push_back(sdf.projectToIsosurface(mesh.V(vertex)-P*u+P*v));
        testx.push_back(sdf.projectToIsosurface(mesh.V(vertex)+P*u+P*v));

    } else { // interior vertex
        // Generate random perturbation points within cube
        for (int i = 0; i < NUM_RAND_PERTURBATIONS; ++i) {
            testx.push_back(mesh.V(vertex)+Vec3f(randf(P),randf(P),randf(P)));
        }
        // and include original point and corners of cube
        testx.push_back(mesh.V(vertex));
        testx.push_back(mesh.V(vertex)+Vec3f(-P,-P,-P));
        testx.push_back(mesh.V(vertex)+Vec3f(+P,-P,-P));
        testx.push_back(mesh.V(vertex)+Vec3f(-P,+P,-P));
        testx.push_back(mesh.V(vertex)+Vec3f(+P,+P,-P));
        testx.push_back(mesh.V(vertex)+Vec3f(-P,-P,+P));
        testx.push_back(mesh.V(vertex)+Vec3f(+P,-P,+P));
        testx.push_back(mesh.V(vertex)+Vec3f(-P,+P,+P));
        testx.push_back(mesh.V(vertex)+Vec3f(+P,+P,+P));
    }

    // Find incident tetrahedra, and the indices that refer to the 
    // current vertex within those tetrahedra.
    std::vector<Tet> incidentTets;
    mesh.getIncidentTets(vertex, incidentTets);
    std::vector<int> incidentVertexIndex(incidentTets.size(), -1);
    for (size_t i = 0; i < incidentTets.size(); ++i) {
        for (int v = 0; v < 4; ++v) {
            if (mesh.V(vertex) == incidentTets[i][v]) {
                incidentVertexIndex[i] = v;
                continue;
            }
        }
        assert(incidentVertexIndex[i] != -1);
    }

	
    // Find quality metrics for tetrahedra created from each perturbation
	// of the given vertex.
	// We're looking for the perturbation that gives the highest quality metric.
	float bestQuality = -FLT_MAX;
	int bestPerturb = -1;

	for (size_t p = 0; p < testx.size(); ++p) {
		// Create tetrahedra to represent the local neighbourhood if 
		// the given vertex were perturbed by the current amount.
		Vec3f perturbedVert = testx[p];
		
		// Compute quality of each tetrahedron and use the worst one
		// as representative of the current perturbation.
		float perturbQuality = FLT_MAX;
		for (size_t t = 0; t < incidentTets.size(); ++t) {

			// Perturb the current vertex in the tetrahedron.
            incidentTets[t][incidentVertexIndex[t]] = perturbedVert;
			
			float tetQuality = compute_tet_quality(incidentTets[t]);
			
			// Keep the lowest-quality tet's quality measure.
			if (tetQuality < perturbQuality) {
				perturbQuality = tetQuality;
			}
		}
		
		// Check whether the current perturbation is the best so far.
		if (perturbQuality > bestQuality) {
			bestQuality = perturbQuality;
			bestPerturb = p;
		} else if (perturbQuality == bestQuality && 
					bestPerturb >= 0) {
			// To break a tie, take the smallest perturbation.
			if (dist2(mesh.V(vertex), testx[p]) < 
                    dist2(mesh.V(vertex), testx[bestPerturb])) {
				bestQuality = perturbQuality;
				bestPerturb = p;
			}
		}
	}
	
	// Pick the perturbation that yields the best quality (lowest value)
	// and apply it to the vertex.
    if (bestPerturb >= 0) {
	    mesh.V(vertex) = testx[bestPerturb];
	    
        // Return the quality metric for the new vertex.
	    return bestQuality;
    }

    // Just in case we failed to get a good perturbation,
    // compute quality of incident tets.
    float quality = FLT_MAX;
    for (size_t t = 0; t < incidentTets.size(); ++t) {
        
        // Move the vertex back to its original location in the tetrahedron.
        incidentTets[t][incidentVertexIndex[t]] = mesh.V(vertex);
        
        float tetQuality = compute_tet_quality(incidentTets[t]);
        
        // Keep the lowest-quality tet's quality measure.
        if (tetQuality < quality) {
            quality = tetQuality;
        }
    }
    return quality;
}


float
optimize_tet_mesh(TetMesh& mesh,
                  const SDF& sdf,
                  const std::vector<int>& boundary_verts,
                  const FeatureSet& featureSet,
                  const std::vector<int>& feature_endpoints,
                  const std::map<int, int>& vertex_feature_map)
{
	std::vector<bool> isBoundary(mesh.vSize(), false);
	std::vector<bool> isFeatureEndpoint(mesh.vSize(), false);
    std::vector<int> vertexFeatures(mesh.vSize(), -1);

    for (size_t b = 0; b < boundary_verts.size(); ++b) {
        assert(boundary_verts[b] < (int)isBoundary.size());
        isBoundary[boundary_verts[b]] = true;
    }
    for (size_t f = 0; f < feature_endpoints.size(); ++f) {
        assert(feature_endpoints[f] < (int)isFeatureEndpoint.size());
        isFeatureEndpoint[feature_endpoints[f]] = true;
    }
    for (std::map<int, int>::const_iterator vfItr = vertex_feature_map.begin();
         vfItr != vertex_feature_map.end(); ++vfItr) {
        assert(vfItr->first < (int)vertexFeatures.size());
        vertexFeatures[vfItr->first] = vfItr->second;
    }

    static const int MAX_ITERATIONS = 20;
    static const float MIN_IMPROVEMENT_RATIO = 1.0;

	// The "quality" of the tetrahedra will be represented by a
	// number that we want to maximize.
	// Therefore, a higher "quality" value => better tets
	float tetQuality = -1.0;

    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        float currTetQuality = FLT_MAX;

        // First loop over interior vertices
        for (size_t v = 0; v < mesh.vSize(); ++v) {
            if (!isBoundary[v] && !isFeatureEndpoint[v]
                && vertexFeatures[v] == -1) {
                float quality = optimize_vertex(v, mesh, sdf, featureSet,
                                                false, -1);
                if (quality < currTetQuality) {
                    currTetQuality = quality;
                }
            }
        }

        // Then loop over boundary vertices
        for (size_t v = 0; v < mesh.vSize(); ++v) {
            if (isBoundary[v] && !isFeatureEndpoint[v] 
                && vertexFeatures[v] == -1) {
                float quality = optimize_vertex(v, mesh, sdf, featureSet,
                                                true, -1);
                if (quality < currTetQuality) {
                    currTetQuality = quality;
                }
            }
        }

        // Then loop over feature vertices.
        for (size_t v = 0; v < mesh.vSize(); ++v) {
            if (vertexFeatures[v] != -1 && !isFeatureEndpoint[v]) {
                float quality = optimize_vertex(v, mesh, sdf, featureSet,
                                                false, vertexFeatures[v]);
                if (quality < currTetQuality) {
                    currTetQuality = quality;
                }
            }
        }

        std::cout << "." << std::flush;

        // Terminate if quality did not improve significantly.
        if (currTetQuality <= tetQuality*MIN_IMPROVEMENT_RATIO) {
            std::cout << std::endl;
            return currTetQuality;
        }

        tetQuality = currTetQuality;
    }
    
    std::cout << std::endl;
    
    return tetQuality;
}


/*
float
optimize_tet_mesh(TetMesh& mesh, 
			      const SDF& sdf,
                  const std::vector<int> &boundary_verts, 
                  const std::vector<int> &feature_verts,
                  const std::map<int, Vec3f> &feature_edge_verts,
                  const std::map<int, int> &vertex_feature_map)
{
	// Do a breadth-first search through the tet mesh, starting at the
	// feature vertices and boundary vertices
	// (these are the ones that have been perturbed out of lattice locations).
	// At each tet vertex, try a bunch of new positions and pick the one
	// that optimizes our tet quality metric for the incident tets.
	
	std::queue<int> bfsQueue;
	std::vector<bool> alreadyTraversed(mesh.v.size(), false);
	std::vector<bool> isBoundary(mesh.v.size(), false);
	std::vector<bool> isFeature(mesh.v.size(), false);
    std::vector<bool> isFeatureEdge(mesh.v.size(), false);
	
	// The "quality" of the tetrahedra will be represented by some
	// number that we want to minimize.
	// Therefore, a lower "quality" value => better tets
	float maxQualityMetric = 0.0f;

    for (std::vector<int>::const_iterator fItr = feature_verts.begin();
            fItr != feature_verts.end(); ++fItr) {
        isFeature[*fItr] = true;
    }
    for (std::map<int, Vec3f>::const_iterator fItr = feature_edge_verts.begin();
            fItr != feature_edge_verts.end(); ++fItr) {
        isFeatureEdge[fItr->first] = true;
    }
    for (std::vector<int>::const_iterator bItr = boundary_verts.begin();
            bItr != boundary_verts.end(); ++bItr) {
        isBoundary[*bItr] = true;
    }
	
	// Add feature vertices to queue, followed by boundary vertices.
	// This ensures that all boundary vertices are visited first by the BFS.
	if (!feature_verts.empty()) {
		for (std::vector<int>::const_iterator fItr = feature_verts.begin();
				fItr != feature_verts.end(); ++fItr) {
            if (!alreadyTraversed[*fItr]) {
			    bfsQueue.push(*fItr);
			    alreadyTraversed[*fItr] = true;
            }
		}
	} else {
        // If there are no feature vertices, start from boundary vertices.
		for (std::vector<int>::const_iterator bItr = boundary_verts.begin();
				bItr != boundary_verts.end(); ++bItr) {
			if (!alreadyTraversed[*bItr]) {
				bfsQueue.push(*bItr);
				alreadyTraversed[*bItr] = true;
			}
		}
	}
	
	// Breadth-First Search loop
	while (!bfsQueue.empty()) {
		int currVertex = bfsQueue.front();
		bfsQueue.pop();
		
		std::vector<int> incidentTets;
		std::set<int> neighbours;
		
		// Find the tets that contain the current vertex,
		// and its neighbouring vertices.
		// Shouldn't need to loop over all tets to get vertex+neighbours...
		for (size_t t = 0; t < mesh.t.size(); ++t) {
			if (mesh.t[t][0] == currVertex || mesh.t[t][1] == currVertex ||
				mesh.t[t][2] == currVertex || mesh.t[t][3] == currVertex) {
				incidentTets.push_back(t);
				neighbours.insert(mesh.t[t][0]);
				neighbours.insert(mesh.t[t][1]);
				neighbours.insert(mesh.t[t][2]);
				neighbours.insert(mesh.t[t][3]);
			}
		}
		
		// Because we didn't guard against the current vertex being added
		// to the neighbours set, remove it from the set.
		neighbours.erase(currVertex);
		
		// Add neighbours to queue if they haven't been already.
		for (std::set<int>::iterator nItr = neighbours.begin();
			 nItr != neighbours.end(); ++nItr) {
			if (!alreadyTraversed[*nItr]) {
				bfsQueue.push(*nItr);
				alreadyTraversed[*nItr] = true;
			}
		}
		
		// Don't move feature vertices.
		if (isFeature[currVertex]) {
			continue;
		}

        // Retrieve feature edge, if it exists.
        Vec3f edge(0.0f, 0.0f, 0.0f);
        std::map<int,Vec3f>::const_iterator e = 
            feature_edge_verts.find(currVertex);
        if (e != feature_edge_verts.end()) {
            edge = e->second;
        }

        // Use vertex_feature_map to know which verts were snapped to
        // which polyline features.
		
		// Do optimization of current vertex location.
		float quality = optimize_vertex(currVertex, mesh, sdf, 
										incidentTets, 
                                        isBoundary[currVertex],
                                        isFeatureEdge[currVertex],
                                        edge);
		
		// Update global quality metric.
		if (quality > maxQualityMetric) {
			maxQualityMetric = quality;
		}
	}
	
	return maxQualityMetric;
}
*/

