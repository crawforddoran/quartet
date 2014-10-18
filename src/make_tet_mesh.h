#ifndef MAKE_TET_MESH_H
#define MAKE_TET_MESH_H

#include "sdf.h"
#include "tet_mesh.h"
#include "feature.h"

// From a signed distance function sdf, create a tetrahedral mesh of 
// its interior.
// We use an isosurface stuffing approach starting from an acute lattice.
//  mesh - The tetrahedral mesh to be constructed.
//  sdf - The signed distance field from which the mesh will be built
//  optimize - if true, the mesh will be optimized after construction
//  outputIntermediateMeshes - if true, the mesh will be output at various 
//                             intermediate stages of the mesh generation 
//                             process.
void
make_tet_mesh(TetMesh& mesh,
              const SDF& sdf,
              bool optimize = true,
              bool outputIntermediateMeshes = false,
              bool unsafeFeatureMatching = false);


// Create tetrahedral mesh as above that also incorporates the 
// sharp features found in the given featureSet.
//  sdf - The signed distance field from which the mesh will be built
//  mesh - The tetrahedral mesh to be constructed.
//  featureSet - A set of features (points and edges) that will be
//               matched as closely as possible by the constructed tet mesh.
//  optimize - if true, the mesh will be optimized after construction
//  outputIntermediateMeshes - if true, the mesh will be output at various 
//                             intermediate stages of the mesh generation 
//                             process.
void 
make_tet_mesh(TetMesh& mesh,
              const SDF& sdf,
              FeatureSet& featureSet,
              bool optimize = true,
              bool outputIntermediateMeshes = false,
              bool unsafeFeatureMatching = false);

#endif
