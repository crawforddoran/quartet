#ifndef SDF_H
#define SDF_H

#include "array3.h"
#include "vec.h"


//
// Signed Distance Field
//
// Scalar field storing the shortest distance to an implicit surface
// (which is therefore found at the zero-isosurface).
// The field has a negative sign inside the surface and positive sign outside.
// The scalar field phi is stored as a 3-dimensional array of values, and
// trilinear interpolation is provided for arbitrary spatial lookups.
//
struct SDF
{
    Array3f phi;
    const Vec3f origin;
    const float dx;

    // Constructor
    // origin is the origin of the 3D array.
    // dx is the grid spacing in 3D space
    // ni*nj*nk are the dimensions of the grid.
    SDF(const Vec3f& origin,
        const float dx,
        int ni, int nj, int nk);

    // Evaluate the scalar field at the given point using 
    // trilinear interpolation.
    float operator() (const Vec3f& x) const;

    // Compute the gradient vector of the field at the given point.
    Vec3f gradient(const Vec3f& x) const;

    // Compute the "normal" vector of the field at the given point.
    // This is equivalent to the normalized gradient vector.
    Vec3f normal(const Vec3f& x) const;

    // Project the given point x onto the zero-isosurface of the field.
    // Note that this is only accurate if x is already near the isosurface.
    Vec3f projectToIsosurface(const Vec3f& x) const;

private:
    const float over_dx; // 1/dx
};

#endif
