#ifndef MATCH_FEATURES_H
#define MATCH_FEATURES_H

#include "tet_mesh.h"
#include "feature.h"
#include "vec.h"
#include <vector>
#include <map>


// Modify the given TetMesh to match the features in the given FeatureSet.
//  mesh - the TetMesh to be modified
//  featureSet - the set of features to be matched by the mesh
//  dx - the resolution of the mesh.  Features below this resolution are not
//       guaranteed to be matched accurately.
//  boundary_verts - indices of the boundary vertices of the mesh
//  boundary_tris - triangles of the boundary mesh, stored as indices into the
//                  vertices of the mesh
//  feature_endpoints - on exit, stores the indices of the mesh vertices that
//                      are now constrained to feature endpoints.
//  vertex_feature_map - on exit, stores a mapping of vertex indices to the
//                       indices of the features to which they are now 
//                       constrained
//  unsafe - if true, allow feature matching operations that might invert
//           tetrahedra
// Returns false iff the endpoints of the features could not be resolved,
// indicating algorithm failure.
bool
match_features(TetMesh& mesh,
               FeatureSet& featureSet,
               float dx,
               const std::vector<int>& boundary_verts,
               const std::vector<Vec3i>& boundary_tris,
               std::vector<int>& feature_endpoints,
               std::map<int, int>& vertex_feature_map,
               bool unsafe = false);

#endif
