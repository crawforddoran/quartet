#ifndef MAKE_SIGNED_DISTANCE_H
#define MAKE_SIGNED_DISTANCE_H

#include "sdf.h"
#include "vec.h"

// Create a signed distance field for a water-tight triangle mesh (absolute
// distances will be pretty good even for arbitrary triangle soup, but
// correct signs rely on having a closed mesh).
// Input tri is a list of triangles in the mesh, and x is the positions of the
// vertices;
// Variable sdf should already be initialized with the desired origin, grid
// size and grid dimensions.  
// Its values will be set up with the signed distance field.
void
make_signed_distance(const std::vector<Vec3i> &tri,
                     const std::vector<Vec3f> &x,
                     SDF &sdf);

#endif
