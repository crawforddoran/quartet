#include "make_tet_mesh.h"
#include "match_features.h"
#include "optimize_tet_mesh.h"
#include "read_obj.h"
#include "geometry_queries.h"
#include "vec.h"
#include <map>
#include <set>
#include <cfloat>
#include <ctime>

// The basic acute tetrahedral lattice definition.
// See Brent Williams' MSc thesis at UBC (supervisor Robert Bridson) for more
// information on where this is from and why it's especially useful for mesh
// generation.
const int num_lattice_nodes=27;
const Vec3i lattice_node[27]={
    Vec3i(1,0,4),
    Vec3i(3,0,4),
    Vec3i(5,0,4),
    Vec3i(0,0,2),
    Vec3i(2,1,2),
    Vec3i(4,0,2),
    Vec3i(1,0,0),
    Vec3i(3,0,0),
    Vec3i(5,0,0),
    Vec3i(0,2,5),
    Vec3i(2,2,4),
    Vec3i(4,2,5),
    Vec3i(0,2,3),
    Vec3i(2,3,2),
    Vec3i(4,2,3),
    Vec3i(0,2,1),
    Vec3i(2,2,0),
    Vec3i(4,2,1),
    Vec3i(1,4,4),
    Vec3i(3,4,4),
    Vec3i(5,4,4),
    Vec3i(0,4,2),
    Vec3i(2,5,2),
    Vec3i(4,4,2),
    Vec3i(1,4,0),
    Vec3i(3,4,0),
    Vec3i(5,4,0),
};
const int num_lattice_tets=46;
const Vec4i lattice_tet[46]={
    Vec4i(12,4,0,3),
    Vec4i(0,12,9,10),
    Vec4i(4,12,0,10),
    Vec4i(1,4,0,10),
    Vec4i(13,12,4,10),

    Vec4i(14,1,2,11),
    Vec4i(1,14,4,10),
    Vec4i(14,1,4,5),
    Vec4i(14,13,4,10),
    Vec4i(1,14,2,5),
    Vec4i(1,11,14,10),

    Vec4i(6,7,16,4),
    Vec4i(6,4,15,3),
    Vec4i(15,4,12,3),
    Vec4i(12,13,4,15),
    Vec4i(16,15,6,4),
    Vec4i(13,15,16,4),

    Vec4i(4,17,7,16),
    Vec4i(17,4,7,5),
    Vec4i(17,13,16,4),
    Vec4i(4,14,5,17),
    Vec4i(8,5,17,7),
    Vec4i(13,17,14,4),

    Vec4i(22,19,18,13),
    Vec4i(18,12,10,9),
    Vec4i(13,19,18,10),
    Vec4i(13,12,18,21),
    Vec4i(12,13,18,10),
    Vec4i(18,22,13,21),
    
    Vec4i(19,13,14,10),
    Vec4i(19,11,10,14),
    Vec4i(11,19,20,14),
    Vec4i(14,19,20,23),
    Vec4i(19,22,23,13),
    Vec4i(13,19,14,23),

    Vec4i(13,12,21,15),
    Vec4i(24,22,21,13),
    Vec4i(15,24,13,16),
    Vec4i(24,15,13,21),
    Vec4i(24,25,22,13),
    Vec4i(25,24,16,13),

    Vec4i(25,22,13,23),
    Vec4i(13,17,16,25),
    Vec4i(25,13,17,23),
    Vec4i(13,14,17,23),
    Vec4i(23,26,17,25)
};


// cut a chunk of the lattice out to roughly match the shape
static void
cut_lattice(const SDF& sdf,
            TetMesh& mesh,
            std::vector<float>& vphi) // phi evaluated at vertices
{
    std::map<Vec3i,int> lattice_map;
    // loop over tiles that overlap bounding box
    int i, j, k;
    for (k=0; k<sdf.phi.nk; k+=4) for (j=0; j<sdf.phi.nj; j+=4) 
      for (i=0; i<sdf.phi.ni; i+=4) {
        for (int t=0; t<num_lattice_tets; ++t) {
            // check this tet is entirely inside the grid & has a negative phi
            bool good_tet = true;
            bool negative_phi = false;
            for (int u=0; u<4; ++u) {
                int a = lattice_node[lattice_tet[t][u]][0]+i;
                if (a >= sdf.phi.ni) { good_tet=false; break; }
                int b = lattice_node[lattice_tet[t][u]][1]+j;
                if (b >= sdf.phi.nj) { good_tet=false; break; }
                int c = lattice_node[lattice_tet[t][u]][2]+k;
                if (c >= sdf.phi.nk) { good_tet=false; break; }
                if (sdf.phi(a,b,c) < 0)
                    negative_phi = true;
            }
            if (good_tet && negative_phi) {
                // add vertices if they don't exist yet, find actual tet to add
                Vec4i actual_tet(-1,-1,-1,-1);
                for (int u=0; u<4; ++u) {
                    Vec3i vertex=lattice_node[lattice_tet[t][u]]+Vec3i(i,j,k);
                    std::map<Vec3i,int>::iterator p = lattice_map.find(vertex);
                    if (p == lattice_map.end()) {
                        actual_tet[u] = mesh.vSize();
                        lattice_map[vertex] = actual_tet[u];
                        mesh.verts().push_back(sdf.origin+sdf.dx*Vec3f(vertex));
                        vphi.push_back(sdf.phi(vertex[0],vertex[1],vertex[2]));
                    } else {
                        actual_tet[u]=p->second;
                    }
                }
                mesh.tets().push_back(actual_tet);
            }
        }
    }
}


// Fix some edge crossings by warping vertices if it's admissible
static void
warp_vertices(const float threshold,
              TetMesh& mesh,
              std::vector<float>& vphi)
{
    assert(threshold>=0 && threshold<=0.5);
    std::vector<float> warp(mesh.vSize(), FLT_MAX);
    std::vector<int> warp_nbr(mesh.vSize(), -1);
    std::vector<Vec3f> d(mesh.vSize(), Vec3f(0,0,0));
    // it's wasteful to iterate through tets just to look at edges; oh well
    for (size_t t=0; t<mesh.tSize(); ++t) {
        for (int u=0; u<3; ++u) {
            int i = mesh.T(t)[u];
            for (int v=u+1; v<4; ++v) {
                int j = mesh.T(t)[v];
                if ((vphi[i]<0 && vphi[j]>0) || (vphi[i]>0 && vphi[j]<0)) {
                    float alpha = vphi[i]/(vphi[i]-vphi[j]);
                    if (alpha < threshold) { // warp i?
                        float d2 = alpha*dist2(mesh.V(i), mesh.V(j));
                        if (d2 < warp[i]) {
                            warp[i] = d2;
                            warp_nbr[i] = j;
                            d[i] = alpha*(mesh.V(j)-mesh.V(i));
                        }
                    } else if (alpha > 1-threshold) { // warp j?
                        float d2 = (1-alpha)*dist2(mesh.V(i), mesh.V(j));
                        if (d2 < warp[j]) {
                            warp[j] = d2;
                            warp_nbr[j] = i;
                            d[j] = (1-alpha)*(mesh.V(i)-mesh.V(j));
                        }
                    }
                }
            }
        }
    }
    // do the warps (also wasteful to loop over all vertices; oh well)
    for (size_t i=0; i<mesh.vSize(); ++i) if(warp_nbr[i]>=0) {
        mesh.V(i) += d[i];
        vphi[i] = 0; 
    }
}


// assuming vphi[i] and vphi[j] have different signs, find or create a new
// vertex on the edge between them where phi interpolates to zero.
static int
cut_edge(int i, int j,
         TetMesh& mesh,
         std::vector<float> &vphi,
         std::map<Vec2i,int> &cut_map)
{
    assert((vphi[i]<0 && vphi[j]>0) || (vphi[i]>0 && vphi[j]<0));
    Vec2i edge(min(i,j), max(i,j));
    std::map<Vec2i,int>::iterator p = cut_map.find(edge);
    if (p == cut_map.end()) { // this edge hasn't been cut yet?
        int v = (int)mesh.vSize();
        float alpha = vphi[i]/(vphi[i]-vphi[j]);
        mesh.verts().push_back((1-alpha)*mesh.V(i) + alpha*mesh.V(j));
        vphi.push_back(0);
        cut_map[edge] = v;
        return v;
    } else { // already done this edge
        return p->second;
    }
}


// Find all tets with vertices where a corner has vphi>0 and trim them back.
static void
trim_spikes(TetMesh& mesh,
            std::vector<float> &vphi)
{
    // keep track of edges we cut - store the index of the new vertex
    std::map<Vec2i,int> cut_map;
    // we have a separate list for the results of trimming tets
    std::vector<Vec4i> new_tets;
    // go through the existing tets to see what we have to trim
    for (int t=0; t<(int)mesh.tSize(); ++t) {
        int p, q, r, s; assign(mesh.T(t), p, q, r, s);
        if (vphi[p]<=0 && vphi[q]<=0 && vphi[r]<=0 && vphi[s]<=0)
            continue; // this tet doesn't stick out

        // We have a plethora of cases to deal with. Let's sort the vertices
        // by their phi values and break ties with index, remembering if we
        // have switched orientation or not, to cut down on the number of
        // cases to handle. Use a sorting network to do the job.
        bool flipped=false;
        if (vphi[p]<vphi[q] || (vphi[p]==vphi[q] && p<q)) {
            std::swap(p, q);
            flipped = !flipped;
        }
        if (vphi[r]<vphi[s] || (vphi[r]==vphi[s] && r<s)) {
            std::swap(r, s);
            flipped = !flipped;
        }
        if (vphi[p]<vphi[r] || (vphi[p]==vphi[r] && p<r)) {
            std::swap(p, r);
            flipped = !flipped;
        }
        if (vphi[q]<vphi[s] || (vphi[q]==vphi[s] && q<s)) {
            std::swap(q, s);
            flipped = !flipped;
        }
        if (vphi[q]<vphi[r] || (vphi[q]==vphi[r] && q<r)) {
            std::swap(q, r);
            flipped = !flipped;
        }

        // sanity checks
        assert(vphi[p]>=vphi[q] && vphi[q]>=vphi[r] && vphi[r]>=vphi[s]); //sort
        assert(vphi[p]>0); // we already skipped interior tets
        assert(vphi[s]<=0); // we never generate tets with all positive phi

        // now do the actual trimming
        if (vphi[s]==0) { // +++0 entirely outside
            mesh.T(t)=mesh.tets().back(); // overwrite with last tet
            mesh.tets().pop_back(); // get rid of the last one
            --t; // decrement to cancel the next for-loop increment

        } else if (vphi[r]>0) { // +++- just one vertex inside, three out
            // replace this tet with one clipped back to the isosurface
            int ps = cut_edge(p, s, mesh, vphi, cut_map),
                qs = cut_edge(q, s, mesh, vphi, cut_map),
                rs = cut_edge(r, s, mesh, vphi, cut_map);
            if (flipped)
                mesh.T(t) = Vec4i(qs, ps, rs, s);
            else
                mesh.T(t) = Vec4i(ps, qs, rs, s);

        } else if(vphi[q]<0) { // +--- just one vertex outside, three in
            // Have to tetrahedralize the resulting triangular prism.
            // Note that the quad faces have to be split in a way that will
            // be consistent with face-adjacent tets: our sort of pqrs with
            // consistently broken ties makes this work. We cut the quad to the
            // deepest vertex (phi as negative as possible).
            int pq = cut_edge(p, q, mesh, vphi, cut_map),
                pr = cut_edge(p, r, mesh, vphi, cut_map),
                ps = cut_edge(p, s, mesh, vphi, cut_map);
            if (flipped) {
                mesh.T(t) = Vec4i(q, pq, r, s);
                new_tets.push_back(Vec4i(r, pq, pr, s));
                new_tets.push_back(Vec4i(s, pq, pr, ps));
            } else {
                mesh.T(t)=Vec4i(pq, q, r, s);
                new_tets.push_back(Vec4i(pq, r, pr, s));
                new_tets.push_back(Vec4i(pq, s, pr, ps));
            }

        } else if (vphi[q]>0 && vphi[r]<0) { // ++-- two vertices out, two in
            int pr = cut_edge(p, r, mesh, vphi, cut_map),
                ps = cut_edge(p, s, mesh, vphi, cut_map),
                qr = cut_edge(q, r, mesh, vphi, cut_map),
                qs = cut_edge(q, s, mesh, vphi, cut_map);
            if (flipped) {
                mesh.T(t) = Vec4i(qr, pr, r, s);
                new_tets.push_back(Vec4i(qs, pr, qr, s));
                new_tets.push_back(Vec4i(ps, pr, qs, s));
            } else {
                mesh.T(t) = Vec4i(pr, qr, r, s);
                new_tets.push_back(Vec4i(pr, qs, qr, s));
                new_tets.push_back(Vec4i(pr, ps, qs, s));
            }

        } else if (vphi[q]==0 && vphi[r]==0) { // +00- 1 out, 1 in, 2 surface
            int ps = cut_edge(p, s, mesh, vphi, cut_map);
            if (flipped)
                mesh.T(t) = Vec4i(q, ps, r, s);
            else
                mesh.T(t) = Vec4i(ps, q, r, s);

        } else if (vphi[q]==0) { // +0-- 1 out, 2 in, 1 surface
            int pr = cut_edge(p, r, mesh, vphi, cut_map),
                ps = cut_edge(p, s, mesh, vphi, cut_map);
            if (flipped) {
                mesh.T(t) = Vec4i(q, pr, r, s);
                new_tets.push_back(Vec4i(q, ps, pr, s));
            } else {
                mesh.T(t) = Vec4i(pr, q, r, s);
                new_tets.push_back(Vec4i(ps, q, pr, s));
            }

        } else { // ++0- two out, one in, one surface
            int ps = cut_edge(p, s, mesh, vphi, cut_map),
                qs = cut_edge(q, s, mesh, vphi, cut_map);
            if (flipped)
                mesh.T(t) = Vec4i(qs, ps, r, s);
            else
                mesh.T(t) = Vec4i(ps, qs, r, s);
        }
    }
    // append all the remaining new tets
    for (size_t t=0; t<new_tets.size(); ++t)
        mesh.tets().push_back(new_tets[t]);
}


// remove tets with corners where phi=0 but are outside the volume
static void
remove_exterior_tets(TetMesh& mesh,
                     const std::vector<float> &vphi,
                     const SDF& sdf)
{
    for (size_t t=0; t<mesh.tSize(); ) {
        int i, j, k, l; assign(mesh.T(t), i, j, k, l);
        assert(vphi[i]<=0 && vphi[j]<=0 && vphi[k]<=0 && vphi[l]<=0);
        if (vphi[i]==0 && vphi[j]==0 && vphi[k]==0 && vphi[l]==0
            && sdf((mesh.V(i)+mesh.V(j)+mesh.V(k)+mesh.V(l))/4)>0) {
            mesh.T(t) = mesh.tets().back();
            mesh.tets().pop_back();
        } else {
            ++t;
        }
    }
}


// form topological structures of use in mesh optimization
/*
void
find_neighborhoods(const int nv, // number of vertices
                   const std::vector<Vec4i>& tet,
                   )
*/                  


float
max_dihedral_angle(const TetMesh& mesh)
{
	// Measure maximum dihedral angle of tet mesh.
	float maxAngle = 0.0f;
	for (size_t t = 0; t < mesh.tSize(); ++t) {
        std::vector<float> angles;
        mesh.getTet(t).dihedralAngles(angles);
        for (size_t i = 0; i < angles.size(); ++i) {
		    if (angles[i] > maxAngle) {
			    maxAngle = angles[i];
		    }
        }
	}
	return maxAngle * 180.0f / M_PI;
}


void
make_tet_mesh(TetMesh& mesh,
              const SDF& sdf,
              bool optimize,
              bool intermediate,
              bool unsafe)
{
    // Create empty placeholder for unused FeatureSet.
    FeatureSet features;

    // Forward function call
    make_tet_mesh(mesh, sdf, features, optimize, intermediate, unsafe);
}


void
make_tet_mesh(TetMesh& mesh,
              const SDF& sdf,
              FeatureSet& featureSet,
              bool optimize,
              bool intermediate,
              bool unsafe)
{
    // Initialize exact arithmetic
    initialize_exact();

    // Initialize mesh from acute lattice.
    std::vector<float> vphi; // SDF value at each lattice vertex.
    mesh.verts().resize(0);
    mesh.tets().resize(0);

    std::cout<<"  cutting from lattice"<<std::endl;
    cut_lattice(sdf, mesh, vphi);

    if (intermediate) {
        static const char* cutLatticeFile = "1_cut_lattice.tet";
        static const char* cutLatticeInfo = "1_cut_lattice.info";
        mesh.writeToFile(cutLatticeFile);
        mesh.writeInfoToFile(cutLatticeInfo);
        std::cout << "Mesh written to " << cutLatticeFile << std::endl;
    }

    // Warp any vertices close enough to phi=0 to fix sign crossings.
    static const float WARP_THRESHOLD = 0.3;
    std::cout<<"  warping vertices"<<std::endl;
    warp_vertices(WARP_THRESHOLD, mesh, vphi);

    if (intermediate) {
        static const char* warpVerticesFile = "2_warp_vertices.tet";
        static const char* warpVerticesInfo = "2_warp_vertices.info";
        mesh.writeToFile(warpVerticesFile);
        mesh.writeInfoToFile(warpVerticesInfo);
        std::cout << "Mesh written to " << warpVerticesFile << std::endl;
    }

    // Cut through tets that still poke out of the level set
    std::cout<<"  trimming spikes"<<std::endl;
    trim_spikes(mesh, vphi);

    if (intermediate) {
        static const char* trimSpikesFile = "3_trim_spikes.tet";
        static const char* trimSpikesInfo = "3_trim_spikes.info";
        mesh.writeToFile(trimSpikesFile);
        mesh.writeInfoToFile(trimSpikesInfo);
        std::cout << "Mesh written to " << trimSpikesFile << std::endl;
    }

    // At this point, there could be some bad tets with all four vertices on
    // the surface but which lie outside the level set.  Get rid of them.
    std::cout<<"  removing exterior tets"<<std::endl;
    remove_exterior_tets(mesh, vphi, sdf);

    // Compact mesh and clear away unused vertices.
    std::cout<<"  compacting mesh"<<std::endl;
    mesh.compactMesh();

    if (intermediate) {
        static const char* removeExtFile = "4_remove_exterior.tet";
        static const char* removeExtInfo = "4_remove_exterior.info";
        static const char* removeExtObj = "4_remove_exterior.obj";
        mesh.writeToFile(removeExtFile);
        mesh.writeInfoToFile(removeExtInfo);
        std::vector<Vec3f> objVerts;
        std::vector<Vec3i> objTris;
        mesh.getBoundary(objVerts, objTris);
        write_objfile(objVerts, objTris, removeExtObj);
        std::cout << "Mesh written to " << removeExtFile << std::endl;
    }

    // Compute maximum dihedral angle of mesh (unoptimized, no features).
    std::cout << "  Maximum dihedral angle = " << max_dihedral_angle(mesh) 
        << std::endl;

    if (featureSet.numFeatures() > 0 || optimize)
    {
        // Identify boundary vertices
        std::vector<int> boundary_verts;
        std::vector<Vec3i> boundary_tris;
        std::cout << "  identifying boundary" << std::endl;
        mesh.getBoundary(boundary_verts, boundary_tris);
        assert(!boundary_verts.empty());
        
        // Snap vertices to given features
        std::vector<int> feature_endpoints;
        std::map<int, int> vertex_feature_map;
        if (featureSet.numFeatures() > 0)
        {
            // Move vertices to match desired features
            std::cout << "  matching features" << std::endl;
            clock_t featureTime = clock();

            match_features(mesh, featureSet, sdf.dx,
                           boundary_verts, boundary_tris,
                           feature_endpoints, vertex_feature_map,
                           unsafe);

            featureTime = clock() - featureTime;
        
            std::cout << "  features matched in " 
                      << ((float)featureTime)/CLOCKS_PER_SEC 
                      << " seconds." << std::endl;
            std::cout << "  Maximum dihedral angle = "
                      << max_dihedral_angle(mesh) << std::endl;

            if (optimize && intermediate)
            {
                static const char* matchFeaturesFile = "5_match_features.tet";
                static const char* matchFeaturesInfo = "5_match_features.info";
                mesh.writeToFile(matchFeaturesFile);
                mesh.writeInfoToFile(matchFeaturesInfo);
                std::cout << "Mesh written to " << matchFeaturesFile 
                          << std::endl;
            }
        } 

        // Finish by optimizing the tetmesh, if desired.
        if (optimize) 
        {
            std::cout << "  optimizing mesh" << std::endl;
            clock_t optTime = clock();

            optimize_tet_mesh(mesh, sdf, 
                              boundary_verts, 
                              featureSet,
                              feature_endpoints,
                              vertex_feature_map);

            optTime = clock() - optTime;
            std::cout << "  Mesh optimization completed in " 
                      << ((float)optTime)/CLOCKS_PER_SEC
                      << " seconds." << std::endl;
            std::cout<< "  Maximum dihedral angle = "
                     << max_dihedral_angle(mesh) << std::endl;
        }

        // DEBUGGING
        // Check for inverted tets
        for (size_t t = 0; t < mesh.tSize(); ++t)
        {
            Tet tet = mesh.getTet(t);
            float signedVolume = tet.volume();
            if (signedVolume <= 0.0)
            {
                std::cerr << "Tet #" << t << " is inverted! " <<
                    "Volume = " << signedVolume << "; " <<
                    "Aspect = " << tet.aspectRatio() << std::endl <<
                    "{" << tet[0] << "} {" << tet[1] << "} {" << tet[2] <<
                    "} {" << tet[3] << "}" << std::endl;
            }
        }
    }
}


