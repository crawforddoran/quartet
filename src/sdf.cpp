#include "sdf.h"


SDF::SDF(const Vec3f &origin_,
         const float dx_,
         int ni, int nj, int nk)
    : origin(origin_), dx(dx_), over_dx(1.0f/dx_)
{
    assert(ni >= 2 && nj >= 2 && nk >= 2);

    // Initialize phi to be upper bound on distance.
    phi.assign(ni, nj, nk, (ni+nj+nk)*dx);
}


float
SDF::operator()(const Vec3f& x) const
{
    float fi = (x[0] - origin[0])*over_dx,
          fj = (x[1] - origin[1])*over_dx,
          fk = (x[2] - origin[2])*over_dx;

    float u = fi - std::floor(fi),
          v = fj - std::floor(fj),
          w = fk - std::floor(fk);

    int i = static_cast<int>(fi), 
        j = static_cast<int>(fj), 
        k = static_cast<int>(fk);

    if (i < 0) { i = 0; u = 0; } else if (i > phi.ni-2) { i = phi.ni-2; u=1; }
    if (j < 0) { j = 0; v = 0; } else if (j > phi.nj-2) { j = phi.nj-2; v=1; }
    if (k < 0) { k = 0; w = 0; } else if (k > phi.nk-2) { k = phi.nk-2; w=1; }

    return (1-u)*( (1-v)*( (1-w)*phi(i,j,k)     + w*phi(i,j,k+1) )
                     + v*( (1-w)*phi(i,j+1,k)   + w*phi(i,j+1,k+1) ) )
             + u*( (1-v)*( (1-w)*phi(i+1,j,k)   + w*phi(i+1,j,k+1) )
                     + v*( (1-w)*phi(i+1,j+1,k) + w*phi(i+1,j+1,k+1) ) );
}


Vec3f
SDF::gradient(const Vec3f& x) const
{
    float fi = (x[0] - origin[0])*over_dx,
          fj = (x[1] - origin[1])*over_dx,
          fk = (x[2] - origin[2])*over_dx;

    float u = fi - std::floor(fi),
          v = fj - std::floor(fj),
          w = fk - std::floor(fk);

    int i = static_cast<int>(fi), 
        j = static_cast<int>(fj), 
        k = static_cast<int>(fk);

    if (i < 0) { i = 0; u = 0; } else if (i > phi.ni-2) { i = phi.ni-2; u=1; }
    if (j < 0) { j = 0; v = 0; } else if (j > phi.nj-2) { j = phi.nj-2; v=1; }
    if (k < 0) { k = 0; w = 0; } else if (k > phi.nk-2) { k = phi.nk-2; w=1; }

    // Compute interpolated finite differences in each dimension.
    Vec3f diff;

    diff[0] = -(1-v) * ( (1-w)*phi(i,j,k)     + w*phi(i,j,k+1) )
                 - v * ( (1-w)*phi(i,j+1,k)   + w*phi(i,j+1,k+1) )
              +(1-v) * ( (1-w)*phi(i+1,j,k)   + w*phi(i+1,j,k+1) )
                 + v * ( (1-w)*phi(i+1,j+1,k) + w*phi(i+1,j+1,k+1) );
    diff[1] = (1-u) *( -((1-w)*phi(i,j,k)     + w*phi(i,j,k+1))
                      + ((1-w)*phi(i,j+1,k)   + w*phi(i,j+1,k+1)) )
                + u *( -((1-w)*phi(i+1,j,k)   + w*phi(i+1,j,k+1)) 
                      + ((1-w)*phi(i+1,j+1,k) + w*phi(i+1,j+1,k+1)) );
    diff[2] = (1-u) * ( (1-v)*(phi(i,j,k+1)     - phi(i,j,k)) 
                        + v * (phi(i,j+1,k+1)   - phi(i,j+1,k)) )
                + u * ( (1-v)*(phi(i+1,j,k+1)   - phi(i+1,j,k))
                        + v * (phi(i+1,j+1,k+1) - phi(i+1,j+1,k)) );

    return diff * over_dx;
}


Vec3f
SDF::normal(const Vec3f& x) const
{
    Vec3f g = gradient(x);
    return g / (mag(g) + 1e-30); // Avoid divide-by-zero
}


Vec3f
SDF::projectToIsosurface(const Vec3f& x) const
{
    float s0 = 0, // zero'th step
          d0 = (*this)(x); // zero'th value
    if (d0 == 0)
        return x; // Already on the isosurface.

    Vec3f g = gradient(x); // Search direction
    float g2 = mag2(g);
    if (g2 == 0)
        return x; // Give up if we're stuck.

    float tol = 1e-3*std::fabs(d0);

    float s = clamp(-d0/g2, -dx, dx); // first step
    float d = (*this)(x+s*g); // first value

    // Search outward until a sign change is found
    for (int steps = 0; steps < 5; ++steps) {
        if (std::fabs(d) <= tol) // converged?
            return x+s*g;
        if (d0 < 0 && d > 0) // found sign change?
            goto refine_projection;
        if (d0 > 0 && d < 0) // found sign change?
            goto refine_projection;
        // otherwise search forward
        d0 = d;
        s0 = s;
        s *= 1.25; // expand outwards a bit
        d = (*this)(x+s*g);
    }
    // If we make it here, we failed to find a sign change...
    if (std::fabs(d0) < std::fabs(d))
        return x;
    else
        return x+s*g;

refine_projection:
    // Secant guess at step size
    s = (d*s0 - d0*s) / (d-d0);
#ifndef NDEBUG
    float d_new = (*this)(x+s*g);
    if (std::fabs(d_new) > 100.0*tol && tol > 1e-5) {
        std::cerr << "projectToIsosurface failed? sdf=" << d_new
                  << ", tol=" << tol << ", ratio=" 
                  << std::fabs(d_new)/tol << std::endl
                  << "x={" << x << "}, gradient={" << g << "}" << std::endl;
    }
#endif
    return x+s*g;
}


