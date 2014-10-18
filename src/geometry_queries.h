#ifndef GEOMETRY_QUERIES_H
#define GEOMETRY_QUERIES_H

#include "vec.h"


// Initializes values needed for exact arithmetic.
// Call this before calling orientation predicates!
double initialize_exact();

// find distance x0 is from segment x1-x2
float
point_segment_distance(const Vec3f &x0,
                       const Vec3f &x1,
                       const Vec3f &x2);

// find distance x0 is from triangle x1-x2-x3
float
point_triangle_distance(const Vec3f &x0,
                        const Vec3f &x1,
                        const Vec3f &x2,
                        const Vec3f &x3);

// Calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2) and return an
// SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
int
orientation(double x1, double y1,
            double x2, double y2,
            double &twice_signed_area);

// Robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3), with SOS
// determination of edge and corner cases.
// If true is returned, the barycentric coordinates are set in a,b,c;
// if false is returned, they could be partially overwritten - do not trust.
bool point_in_triangle_2d(double x0, double y0, 
                          double x1, double y1,
                          double x2, double y2,
                          double x3, double y3,
                          double& a, double& b, double& c);


// Calculate six times the signed area of the tetrahedron
// (a_x,a_y,a_z)-(b_x,b_y,b_z)-(c_x,c_y,c_z)-(0,0,0)
// and return an SOS-determined sign
// (-1, +1, or 0 only if the tet is truly degenerate).
int orientation3D(double a_x, double a_y, double a_z,
                  double b_x, double b_y, double b_z,
                  double c_x, double c_y, double c_z,
                  double d_x, double d_y, double d_z,
                  double &six_signed_volume);

#endif
