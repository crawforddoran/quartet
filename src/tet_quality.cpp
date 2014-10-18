#include "tet_quality.h"


static const float SQRT2 = sqrt(2.0);
static const float PI_3_INV = 3.0/M_PI;
static const float OPT_DIHEDRAL = asin(2.0*sqrt(2.0)/3.0);


float
max_dihedral_angle(const Tet& tet)
{
    int orient = tet.orientation();

    if (orient == 0)
    {
        return 0.0;
    }

    std::vector<float> angles;
    tet.dihedralAngles(angles);
    float maxAngle = angles[0];
    for (size_t i = 1; i < angles.size(); ++i)
    {
        if (angles[i] > maxAngle)
        {
            maxAngle = angles[i];
        }
    }

    // maxAngle will be between OPT_DIHEDRAL and M_PI.
    // Normalize to the range (0,1) and flip so that 0 corresponds to M_PI
    // and 1 corresponds to OPT_DIHEDRAL.
    float quality = 1.0 - (maxAngle-OPT_DIHEDRAL)/(M_PI-OPT_DIHEDRAL);

    if (orient < 0)
    {
        quality *= -1.0;
    }

	return quality;
}


float
min_sine_dihedral_angle(const Tet& tet)
{
    int orient = tet.orientation();
    if (orient == 0)
    {
        return 0.0;
    }

    std::vector<float> angles;
    tet.dihedralAngles(angles);
    float minSin = sin(angles[0]);
    for (size_t i = 1; i < angles.size(); ++i)
    {
        float sinAngle = sin(angles[i]);
        if (minSin > sinAngle)
        {
            minSin = sinAngle;
        }
    }

    if (orient < 0)
    {
        minSin *= -1.0;
    }

    return minSin;
}


float
min_face_angle(const Tet& tet)
{
    int orient = tet.orientation();
    if (orient == 0)
    {
        return 0.0;
    }

    std::vector<float> angles;
    tet.faceAngles(angles);
    float minAngle = angles[0];
    for (size_t i = 1; i < angles.size(); ++i)
    {
        if (angles[i] < minAngle)
        {
            minAngle = angles[i];
        }
    }

    // Normalize
    minAngle *= PI_3_INV;

    // Compute orientation to determine sign of result.
    if (orient < 0)
    {
        minAngle *= -1.0;
    }

    return minAngle;
}


float
aspect_ratio(const Tet& tet)
{
    // Want to minimize aspect ratio, so invert.
    // Normalize by sqrt(2).
    return SQRT2 / tet.aspectRatio();
}


