#ifndef TET_QUALITY_H
#define TET_QUALITY_H

#include "tet_mesh.h"

//
// Implementations of different quality metric for tetrahedra.
// All quality metrics should conform to the following:
//  - Equilateral tetrahedra should have a quality measure of 1.0, which is 
//    the maximum (ie: quality measure should be normalized).
//  - Degenerate tetrahedra should have a quality measure of zero
//  - Inverted tetrahedra should have a negative quality measure.
//
// Ideally these measures should also be scale-invariant.
//
// This gets its own file (ie: is not part of the Tet or TetMesh classes)
// because we want to be able to use these quality metrics interchangeably.
// One user might want to optimize based on max dihedral angle, whereas
// another might want to use aspect ratio.
//
// So ideally, client code should be able to choose which quality function
// (or maybe function object/functor in the future) they want to use.
// Then any code that calls compute_tet_quality will use the chosen 
// quality function without needing to know about it.
//

//
// Quality metrics
//

// Largest dihedral angle
float
max_dihedral_angle(const Tet& tet);

// Minimum sine of dihedral angle
float
min_sine_dihedral_angle(const Tet& tet);

// Smallest face angle
float
min_face_angle(const Tet& tet);

// Aspect ratio is the longest edge length divided by the shortest altitude
float
aspect_ratio(const Tet& tet);


//
// Abstract quality metric function.
//

// Compute a quality metric for the given tetrahedron.
// The quality metric has an optimal (normalized) value of 1.0.
// Degenerate tets have a quality value of 0.0, and inverted tets have 
// a negative quality value.
inline float
compute_tet_quality(const Tet& tet) {
    return aspect_ratio(tet);

    // TODO: Appropriate design pattern to allow plugging in 
    // different metrics as desired.
}

#endif
