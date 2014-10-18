#ifndef OPTIMIZE_TET_MESH_H
#define OPTIMIZE_TET_MESH_H

#include "tet_mesh.h"
#include "sdf.h"
#include "feature.h"
#include <vector>
#include <map>


// Optimize tetrahedral mesh to improve tet quality.
//  mesh - the tetrahedral mesh to optimize
//  sdf - the signed distance field representing the object the tet mesh
//        should reproduce
//  boundary_verts - the indices of the vertices on the boundary of the
//                   tet mesh
//  featureSet - The set of features that the tet mesh is constrained to match
//  feature_endpoints - the indices of the vertices constrained to feature
//                      endpoints. These vertices will not be moved.
//  feature_edge_verts - maps vertex indices to the feature to which they  
//                       they should be constrained
float
optimize_tet_mesh(TetMesh& mesh,
                  const SDF& sdf,
                  const std::vector<int>& boundary_verts,
                  const FeatureSet& featureSet,
                  const std::vector<int>& feature_endpoints,
                  const std::map<int, int>& vertex_feature_map);

#endif
