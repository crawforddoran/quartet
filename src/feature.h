#ifndef FEATURE_H
#define FEATURE_H

#include "trimesh.h"
#include "vec.h"
#include <vector>
#include <map>

//
// A Feature is a point or connected sequence of edges that we wish to 
// preserve when meshing a 3D model.
//
struct Feature 
{
    // The vertices that define the feature when connected by line segments.
    std::vector<Vec3f> vertices;

    // Returns true iff this feature is a single point.
    bool isPoint() const { return vertices.size() == 1; }

    // Return the endpoints of the feature.
    Vec3f getEndpoint1() const { return vertices[0]; }
    Vec3f getEndpoint2() const { return vertices[vertices.size()-1]; }

    // Returns the point on the feature closest to the given point p.
    Vec3f projectToFeature(const Vec3f &p) const;

    // Project the given point p onto the feature and return the tangent 
    // to the feature at that point.
    Vec3f getTangent(const Vec3f &p) const;
};


//
// FeatureSet stores and maintains a set of features.
//
// The following invariant is enforced at all times:
//   - No two features share the same vertex on their interiors.
// In other words, features may only share vertices if those vertices are 
// their endpoints.
//
// The FeatureSet can consolidate its features by combining vertices that
// are within a given distance threshold.
//
// The FeatureSet can be initialized by passing in Features and lists of 
// Features, or it can be initialized by reading a .FEAT file.
//
// Features can be retrieved from the set by numerical index.
// Note that modifying the FeatureSet in any way can invalidate previous
// indices into the FeatureSet.
//
// In addition, the FeatureSet can provide the set of all endpoints of its 
// features.  These vertices are returned with no duplicates.
//
class FeatureSet
{
public:
    // Default constructor
    FeatureSet() {}

    // Constructor
    // f is a list of Features to add to the set.
    FeatureSet(const std::vector<Feature>& f);

    // Add the given list of features to the set.
    void addFeatures(const std::vector<Feature>& features);

    // Add the given feature to the set.
    void addFeature(const Feature& feature);

    // Initialize this FeatureSet based on the given file.
    // The file should be of the .FEAT format,
    // ie: it should contain lists of vertices each prefaced by an
    // integer indicating the number of vertices in the current feature.
    // Returns true iff the file was read successfully.
    bool readFromFile(const char *filename_format, ...);

    // Write this FeatureSet to a file with the given name.
    // Ideally, this should be a .FEAT file.
    // The output format will match the description above.
    bool writeToFile(const char *filename_format, ...);

    // Combine vertices that are within the given distance threshold apart.
    // If this causes any features to share vertices, split those features
    // at those vertices.
    void consolidate(float dist = 0.0f);

    // Returns the number of features in the set.
    unsigned int numFeatures() const { return indices.size(); }

    // Copies the feature at index i in to the given feature parameter.
    void getFeature(unsigned int i, Feature& feature) const;

    // Copies the features in the set into the given features parameter.
    void getFeatures(std::vector<Feature>& features) const;

    // Copies a list of feature endpoints into the given endpoints parameter.
    // The endpoints will not be duplicated; if a vertex is the endpoint of
    // multiple features it will only be returned once.
    void getFeatureEndpoints(std::vector<Vec3f>& endpoints) const;

    // Merges the two given feature endpoints by moving endpoint2 to endpoint1
    // and replacing all references to endpoint2 with references to endpoint1.
    // The given endpoints must be valid endpoints in the FeatureSet. 
    // Ideally they should come from a call to getFeatureEndpoints, 
    // and should not be modified before being passed to this function.
    void mergeEndpoints(const Vec3f& endpoint1, const Vec3f& endpoint2,
                        std::vector<Vec3f>& newEndpoints);

    // Process the given triangle mesh to find sharp features and add them to
    // this FeatureSet.
    // Uses the given angle_threshold to determine which edges and vertices are
    // "sharp".
    void autoDetectFeatures(const TriMesh& mesh, 
                            float angle_threshold_degrees);

private:
    // The set of vertices in the features of the FeatureSet.
    std::vector<Vec3f> vertices;

    // List of features as indices into the vertex list (vertices).
    // indices[f][v] gives the index of the v'th vertex of the f'th feature.
    std::vector< std::vector<int> > indices;

    // Map of vertex positions to their index in vertices.
    std::map<Vec3f, int> vertexMap;

    // Examine the feature at the given index and check if it forms a loop.
    // A loop is a feature that starts and ends at the same vertex.
    // If the feature is degenerate (ie all vertices are the same), 
    // convert it into a point.
    // If the is a loop and the entire feature lies within the given
    // collapseRadius, collapse the feature to a point.
    // If the feature is a loop and extends outside the collapseRadius,
    // split it in half.
    // Returns true iff the feature is split and a new feature is added.
    bool checkForLoop(int fIndex, float collapseRadiusSquared = 0.0f);

    // Find the features in the set containing the given index in its interior.
    // Split that feature in two at that index.
    // Repeat until the given index is only found at the endpoints of features.
    // This may increase the number of features in the FeatureSet, 
    // but the number of vertices remains the same.
    void splitFeatures(int vIndex);
};

#endif
