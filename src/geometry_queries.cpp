#include <cfloat>
#include "geometry_queries.h"

// External functions defined in predicates.cpp
extern double exactinit();
extern double orient3d(double* pa, double* pb, double* pc, double* pd);


double
initialize_exact()
{
    return exactinit();
}

float
point_segment_distance(const Vec3f &x0,
                       const Vec3f &x1,
                       const Vec3f &x2)
{
    Vec3f dx(x2-x1);
    double m2=mag2(dx);
    // find parameter value of closest point on segment
    float s12=(float)(dot(x2-x0, dx)/m2);
    if(s12<0){
        s12=0;
    }else if(s12>1){
        s12=1;
    }
    // and find the distance
    return dist(x0, s12*x1+(1-s12)*x2);
}

float
point_triangle_distance(const Vec3f &x0,
                        const Vec3f &x1,
                        const Vec3f &x2,
                        const Vec3f &x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    Vec3f x13(x1-x3), x23(x2-x3), x03(x0-x3);
    float m13=mag2(x13), m23=mag2(x23), d=dot(x13,x23);
    float invdet=1.f/max(m13*m23-d*d,1e-30f);
    float a=dot(x13,x03), b=dot(x23,x03);
    // the barycentric coordinates themselves
    float w23=invdet*(m23*a-d*b);
    float w31=invdet*(m13*b-d*a);
    float w12=1-w23-w31;
    if(w23>=0 && w31>=0 && w12>=0){ // if we're inside the triangle
        return dist(x0, w23*x1+w31*x2+w12*x3); 
    }else{ // we have to clamp to one of the edges
        if(w23>0) // this rules out edge 2-3 for us
            return min(point_segment_distance(x0,x1,x2),
                       point_segment_distance(x0,x1,x3));
        else if(w31>0) // this rules out edge 1-3
            return min(point_segment_distance(x0,x1,x2),
                       point_segment_distance(x0,x2,x3));
        else // w12 must be >0, ruling out edge 1-2
            return min(point_segment_distance(x0,x1,x3),
                       point_segment_distance(x0,x2,x3));
    }
}

// Find x+y=a+b, but with |x|>|y| (or both zero) non-overlapping in bits.
static void
two_sum(double a, double b, double& x, double& y)
{
    x=a+b;
    double z=x-a;
    y=(a-(x-z))+(b-z);
}

// Find x+y=a with x and y using only half the mantissa bits.
static void
split(double a, double& x, double& y)
{
    double c=134217729*a;
    x=c-(c-a);
    y=a-x;
}

// Find x+y=a*b (exactly) with |x|>|y| (or both zero) nonoverlapping in bits.
static void
two_product(double a,
            double b,
            double& x,
            double& y)
{
    x=a*b;
    double a1, a2, b1, b2;
    split(a, a1, a2);
    split(b, b1, b2);
    y=a2*b2-(((x-a1*b1)-a2*b1)-a1*b2);
}

int
orientation(double x1, double y1,
            double x2, double y2,
            double &twice_signed_area)
{
    twice_signed_area=y1*x2-x1*y2;
    double bound=5*DBL_EPSILON*(std::fabs(y1*x2)+std::fabs(x1*y2));
    if(std::fabs(twice_signed_area)<bound){ // rounding error alert!
        // We'll do the exact version instead with floating point expansions.
        // This could be faster, but not a critical thing to optimize.
        // First find a+b=y1*x2 and c+d=-x1*y2
        double a, b, c, d;
        two_product(y1, x2, a, b);
        two_product(-x1, y2, c, d);
        // Then add c+(a+b) to get f+g+h
        double e, f, g, h;
        two_sum(c, b, e, h);
        two_sum(e, a, f, g);
        // And finally add d+(f+g+h) to get i+j+l+n
        double i, j, k, l, m, n;
        two_sum(d, h, m, n);
        two_sum(m, g, k, l);
        two_sum(k, f, i, j);
        twice_signed_area=i+(j+(l+n));
    }
    if(twice_signed_area>0) return 1;
    else if(twice_signed_area<0) return -1;
    else if(y2>y1) return 1;
    else if(y2<y1) return -1;
    else if(x1>x2) return 1;
    else if(x1<x2) return -1;
    else return 0; // only true when x1==x2 and y1==y2
}

bool point_in_triangle_2d(double x0, double y0, 
                          double x1, double y1,
                          double x2, double y2,
                          double x3, double y3,
                          double& a, double& b, double& c)
{
    x1-=x0; x2-=x0; x3-=x0;
    y1-=y0; y2-=y0; y3-=y0;
    int signa=orientation(x2, y2, x3, y3, a);
    if(signa==0) return false;
    int signb=orientation(x3, y3, x1, y1, b);
    if(signb!=signa) return false;
    int signc=orientation(x1, y1, x2, y2, c);
    if(signc!=signa) return false;
    double sum=a+b+c;
    assert(sum!=0);
    a/=sum;
    b/=sum;
    c/=sum;
    return true;
}


int orientation3D(double a_x, double a_y, double a_z,
                  double b_x, double b_y, double b_z,
                  double c_x, double c_y, double c_z,
                  double d_x, double d_y, double d_z,
                  double &six_signed_volume)
{
    double a[3], b[3], c[3], d[3];
    a[0] = a_x; a[1] = a_y; a[2] = a_z;
    b[0] = b_x; b[1] = b_y; b[2] = b_z;
    c[0] = c_x; c[1] = c_y; c[2] = c_z;
    d[0] = d_x; d[1] = d_y; d[2] = d_z;
    six_signed_volume = orient3d(a, b, c, d);

    // For now our tets have the opposite orientation from Shewchuk's,
    // so we'll flip the sign of the volume.
    six_signed_volume *= -1.0;

    if (six_signed_volume < 0.0) return -1;
    else if (six_signed_volume > 0.0) return 1;
    else return 0;
}


