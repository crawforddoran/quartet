#include <cstdio>
#include <cstring>
#include <fstream>
#include <set>
#include "gluvi.h"
#include "vec.h"
#include "feature.h"
#include "read_obj.h"
#include "tet_mesh.h"
#include "tet_quality.h"

std::vector<Vec3f> x;
std::vector<Vec4i> tet;
std::set<Vec3i> boundary_face_set;

std::vector<Vec3i> tri;

bool drawTriMesh;

bool cutX = false, cutY = false, cutZ = false;
float px, py, pz;

int worst_n = 0;
std::vector<float> tet_quality;

FeatureSet features;

Vec3f centroid;

// decides if a vertex is cut off by the cutting planes:
// x-px = 0, y-py = 0, z-py = 0
bool isCutOff(const Vec3f& v){
    bool isCut = false;
    if(cutX){
	    if(v[0] > px) isCut = true;
    }
    if(cutY){
	    if(v[1] > py) isCut = true;
    }
    if(cutZ){
	    if(v[2] > pz) isCut = true;
    }
    return isCut;
}

// Decides if the first tet has worse quality than the second
bool isWorse(const Vec4i& tet1, const Vec4i& tet2){
    Tet t1, t2;
    t1.v[0] = x[tet1[0]];
    t1.v[1] = x[tet1[1]];
    t1.v[2] = x[tet1[2]];
    t1.v[3] = x[tet1[3]];
    t2.v[0] = x[tet2[0]];
    t2.v[1] = x[tet2[1]];
    t2.v[2] = x[tet2[2]];
    t2.v[3] = x[tet2[3]];
    return compute_tet_quality(t1) < compute_tet_quality(t2);

}

// Rearrange tets so that the worst N tets are at the front
void sortTets(){
    std::nth_element(tet.begin(), tet.begin() + worst_n -1, tet.end(), isWorse);
}

bool isBoundaryFace(int v1, int v2, int v3)
{
    // Try all permutations.
    Vec3i perms[] = {
        Vec3i(v1, v2, v3),
        Vec3i(v2, v3, v1),
        Vec3i(v3, v1, v2),
        Vec3i(v1, v3, v2),
        Vec3i(v2, v1, v3),
        Vec3i(v3, v2, v1)
    };

    for (int i = 0; i < 6; ++i)
    {
        if (boundary_face_set.find(perms[i]) != boundary_face_set.end())
        {
            return true;
        }
    }

    return false;
}


void display(void)
{
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(0.5, 0.5);

    glTranslatef(-centroid[0], -centroid[1], -centroid[2]);

    // Tet Faces
    if (!drawTriMesh) {
        glBegin(GL_TRIANGLES);
        
        // Draw worst n tets
        if(worst_n > 0){
	        for(int t=0; t< worst_n; ++t){
                Vec3f p=x[tet[t][0]],
	                  q=x[tet[t][1]],
	                  r=x[tet[t][2]],
	                  s=x[tet[t][3]];
	
	            // if any vertex is cut off, don't draw tet
	            if(isCutOff(p) || isCutOff(q) || isCutOff(r) || isCutOff(s))
	            	continue;
	            
                /*
	            Vec3f c=(p+q+r+s)/4;
	            p=0.8*p+0.2*c; 
	            q=0.8*q+0.2*c;
	            r=0.8*r+0.2*c;
	            s=0.8*s+0.2*c;
	            //float h=0.25*(t+0.5)/tet.size();
	            float h = 0.25*(t+0.5)/(float)tet.size();
	            if (t%4 == 0)       glColor3f(h+0.7,h,h);
	            else if (t%4 == 1)  glColor3f(h+0.2,h+0.5,h);
	            else if (t%4 == 2)  glColor3f(h,h+0.5,h+0.2);
	            else                glColor3f(h,h,h+0.7);
                */
                glColor3f(1.0, 0.2, 0.0);
                glNormal3fv(normalized(cross((q-p),(r-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(q.v);
	            glVertex3fv(r.v);
                glNormal3fv(normalized(cross((s-p),(q-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(s.v);
	            glVertex3fv(q.v);
                glNormal3fv(normalized(cross((r-p),(s-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(r.v);
	            glVertex3fv(s.v);
                glNormal3fv(normalized(cross((r-q),(s-q))).v);
	            glVertex3fv(q.v);
	            glVertex3fv(r.v);
	            glVertex3fv(s.v);
	        }
	        
     	} else {
	     	// Draw all tets
	        for(size_t t=0; t<tet.size(); ++t){
	            Vec3f p=x[tet[t][0]],
	                  q=x[tet[t][1]],
	                  r=x[tet[t][2]],
	                  s=x[tet[t][3]];
	
	            // if any vertex is cut off, don't draw tet
	            if(isCutOff(p) || isCutOff(q) || isCutOff(r) || isCutOff(s))
	            	continue;
	            
                /* 
	            Vec3f c=(p+q+r+s)/4;
	            p=0.8*p+0.2*c; 
	            q=0.8*q+0.2*c;
	            r=0.8*r+0.2*c;
	            s=0.8*s+0.2*c;
	            //float h=0.25*(t+0.5)/tet.size();
	            float h = 0.25*(t+0.5)/(float)tet.size();
                */
                /*
	            if (t%4 == 0)       glColor3f(h+0.7,h,h);
	            else if (t%4 == 1)  glColor3f(h+0.2,h+0.5,h);
	            else if (t%4 == 2)  glColor3f(h,h+0.5,h+0.2);
	            else                glColor3f(h,h,h+0.7);
                */

                // Determine colors by whether or not face is on boundary.
                if (isBoundaryFace(tet[t][0], tet[t][1], tet[t][2]))
                {
                    glColor3f(1.0, 1.0, 0.7);
                }
                else
                {
                    glColor3f(0.7, 0.7, 1.0);
                }
                glNormal3fv(normalized(cross((q-p),(r-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(q.v);
	            glVertex3fv(r.v);
                if (isBoundaryFace(tet[t][0], tet[t][3], tet[t][1]))
                {
                    glColor3f(1.0, 1.0, 0.7);
                }
                else
                {
                    glColor3f(0.7, 0.7, 1.0);
                }
                glNormal3fv(normalized(cross((s-p),(q-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(s.v);
	            glVertex3fv(q.v);
                if (isBoundaryFace(tet[t][0], tet[t][2], tet[t][3]))
                {
                    glColor3f(1.0, 1.0, 0.7);
                }
                else
                {
                    glColor3f(0.7, 0.7, 1.0);
                }
                glNormal3fv(normalized(cross((r-p),(s-p))).v);
	            glVertex3fv(p.v);
	            glVertex3fv(r.v);
	            glVertex3fv(s.v);
                if (isBoundaryFace(tet[t][1], tet[t][2], tet[t][3]))
                {
                    glColor3f(1.0, 1.0, 0.7);
                }
                else
                {
                    glColor3f(0.7, 0.7, 1.0);
                }
                glNormal3fv(normalized(cross((r-q),(s-q))).v);
	            glVertex3fv(q.v);
	            glVertex3fv(r.v);
	            glVertex3fv(s.v);
	        }
    	}
        glEnd();

        // Tet Edges
        glDisable(GL_LIGHTING);
        glLineWidth(1.0);
        //glColor3f(1.0, 1.0, 0.0);
        glColor3f(0.2, 0.2, 0.2);
        glBegin(GL_LINES);
        for (size_t t=0; t<tet.size(); ++t) {
            Vec3f p=x[tet[t][0]],
                  q=x[tet[t][1]],
                  r=x[tet[t][2]],
                  s=x[tet[t][3]];
			
            // if any vertex is cut off, don't draw tet
            if(isCutOff(p) || isCutOff(q) || isCutOff(r) || isCutOff(s))
            	continue;
            
            /*
            Vec3f c=(p+q+r+s)/4;
            //p=0.8*p+0.2*c; 
            //q=0.8*q+0.2*c;
            //r=0.8*r+0.2*c;
            //s=0.8*s+0.2*c;
            p=0.81*p+0.19*c;
            q=0.81*q+0.19*c;
            r=0.81*r+0.19*c;
            s=0.81*s+0.19*c;
            */
            
            glVertex3fv(p.v);
            glVertex3fv(q.v);

            glVertex3fv(p.v);
            glVertex3fv(r.v);

            glVertex3fv(p.v);
            glVertex3fv(s.v);

            glVertex3fv(q.v);
            glVertex3fv(r.v);

            glVertex3fv(q.v);
            glVertex3fv(s.v);

            glVertex3fv(r.v);
            glVertex3fv(s.v);
        }
        glEnd();

    } else {
        // Draw triangle mesh
        glColor3f(0.9, 0.9, 0.6);
        glBegin(GL_TRIANGLES);
        for (size_t t = 0; t < tri.size(); ++t) {
            glNormal3fv(normalized(cross(x[tri[t][1]]-x[tri[t][0]],
                                         x[tri[t][2]]-x[tri[t][0]])).v);
            glVertex3fv(x[tri[t][0]].v);
            glVertex3fv(x[tri[t][1]].v);
            glVertex3fv(x[tri[t][2]].v);
        }
        glEnd();

        /*
        // Draw triangle edges
        glDisable(GL_LIGHTING);
        glLineWidth(1.0);
        glColor3f(0.2, 0.2, 0.2);
        glBegin(GL_LINES);
        for (size_t t = 0; t < tri.size(); ++t) {
            glVertex3fv(x[tri[t][0]].v);
            glVertex3fv(x[tri[t][1]].v);
            
            glVertex3fv(x[tri[t][1]].v);
            glVertex3fv(x[tri[t][2]].v);
            
            glVertex3fv(x[tri[t][2]].v);
            glVertex3fv(x[tri[t][0]].v);
        }
        glEnd();
        */
    }

    // Features
    glDisable(GL_LIGHTING);
    glLineWidth(5.0);
    glColor3f(0.9, 0.1, 0.7);
    for (size_t f = 0; f < features.numFeatures(); ++f) {
        Feature feat;
        features.getFeature(f, feat);
        glBegin(GL_LINE_STRIP);
        for (size_t v = 0; v < feat.vertices.size(); ++v) {
            glVertex3fv(feat.vertices[v].v);
        }
        glEnd();
    }

    // Feature Endpoints
    glPointSize(10.0);
    glColor3f(0.0, 0.0, 0.0);
    for (size_t f = 0; f < features.numFeatures(); ++f) {
        Feature feat;
        features.getFeature(f, feat);
        glBegin(GL_POINTS);
        glVertex3fv(feat.getEndpoint1().v);
        if (!feat.isPoint()) {
            glVertex3fv(feat.getEndpoint2().v);
        }
        glEnd();
    }
}

int main(int argc, char **argv)
{
    bool usage = false, failure = false;
    drawTriMesh = false;
    bool feature = false;
    char* featureFile = NULL;

    if(argc<2 || argc>7) {
        usage = true;
        failure = true;
    }

    // Parse command line arguments
    for (int arg = 2; arg < argc && !failure; ++arg)
    {
        if (strcmp(argv[arg], "-s") == 0) {
            drawTriMesh = true;

        } else if (strcmp(argv[arg], "-f") == 0) {
            if (feature) {
                std::printf("ERROR: Feature file already specified.");
                usage = true;
                failure = true;
                break;
            }
            ++arg;
            if (arg >= argc) {
                std::printf("ERROR: No feature file provided after "
                            "-f option\n");
                usage = true;
		    } else if (strcmp(argv[arg], "-p") == 0) {
			std::printf("ERROR: No feature file provided after "
	                            "-f option\n");
                usage = true;
            } else {
                feature = true;
                featureFile = argv[arg];
            }

		} else if(strcmp(argv[arg], "-p") == 0) {
			++arg;
		    if (arg >= argc) {
                std::printf("ERROR: No cutting plane specified after "
                            "-p option\n");
                usage = true;

            } else {
	            
                char cutPlane;
                float cutVal;
                
                if (sscanf(argv[arg], "%c%f", &cutPlane, &cutVal) == 2) {
	                
	                if(cutPlane == 'x'){
		                cutX = true;
		                px = cutVal;
	                } else if(cutPlane == 'y'){
		                cutY = true;
		                py = cutVal;
	                } else if(cutPlane == 'z'){
		                cutZ = true;
		                pz = cutVal;
	                } else{
		                std::printf("ERROR: Wrong format in specifying cutting plane\n");
                		usage = true;
	                }
	                
                } else {
	                std::printf("ERROR: Wrong format in specifying cutting plane\n");
                	usage = true;
	                
                }
            }
            
		} else if(strcmp(argv[arg], "-w") == 0) {
			worst_n = atoi(argv[++arg]);
			if ( worst_n == 0) {
				std::printf("ERROR: A proper integer cannot be found after the '-w' option\n");
                	usage = true;
			}
			
			
        } else { // Invalid argument
            usage = true;
            failure = true;
        }
    }

    if (usage) {
        std::printf("Usage: view_tet <input.tet> [-s] [-f <feature.feat>] [-p (x|y|z)<cutting-plane-value>]* [-w num]\n");
    }
    if (failure) {
        return 1;
    }

    Gluvi::init(argv[0], &argc, argv);
    Gluvi::Target3D cam;
    Gluvi::camera=&cam;


    if (!drawTriMesh) {
        // Read .tet file
        FILE *in=std::fopen(argv[1], "r");
        if (!in) {
            std::printf("Error: could not read tet mesh file: %s\n", argv[1]);
            return 1;
        }

        if(std::fgetc(in)!='t') return 1;
        if(std::fgetc(in)!='e') return 1;
        if(std::fgetc(in)!='t') return 1;
        if(std::fgetc(in)!=' ') return 1;
        int nv=-1, nt=-1;
        std::fscanf(in, "%d %d", &nv, &nt);
        if(nv<=0 || nt<=0) return 2;
        x.resize(nv);
        tet.resize(nt);
        for(int i=0; i<nv; ++i)
            std::fscanf(in, "%g %g %g", &x[i][0], &x[i][1], &x[i][2]);
        for(int i=0; i<nt; ++i)
            std::fscanf(in, "%d %d %d %d", &tet[i][0], &tet[i][1], &tet[i][2], &tet[i][3]);
        std::fclose(in);

        // Compute boundary mesh and populate boundary_face_set
        TetMesh tetMesh(x, tet);
        std::vector<int> boundaryVerts;
        std::vector<Vec3i> boundaryTris;
        tetMesh.getBoundary(boundaryVerts, boundaryTris);
        for (size_t i = 0; i < boundaryTris.size(); ++i)
        {
            boundary_face_set.insert(boundaryTris[i]);
        }
        

        if (!x.size()) {
            std::printf("Error reading tet mesh file %s: "
                        "no vertex data found.", argv[1]);
            return 1;
        }
        
        // Reorder to have worst N tets in front
        if(worst_n > 0){
	        sortTets();
        }
        
    } else {
        bool result = read_objfile(x, tri, argv[1]);
        if (!result) {
            std::printf("Error: could not read tri mesh file: %s\n", argv[1]);
            return 1;
        }
        if (!x.size()) {
            std::printf("Error reading tri mesh file %s: "
                        "no vertex data found.", argv[1]);
            return 1;
        }
    }

    if (feature) {
        // Read feature file
        features.readFromFile(featureFile);
    }

    // Construct bounding box and use to find centroid.
    Vec3f bboxMin = x[0], bboxMax = x[0];
    for (size_t i = 1; i < x.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            if (x[i][j] < bboxMin[j]) bboxMin[j] = x[i][j];
            if (x[i][j] > bboxMax[j]) bboxMax[j] = x[i][j];
        }
    }
    centroid = bboxMin + (bboxMax - bboxMin)*0.5;

    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    Gluvi::userDisplayFunc=display;

    /*
    for(int t=0; t<tet.size(); ++t){
        Vec3f p=x[tet[t][0]];
        Vec3f q=x[tet[t][1]];
        Vec3f r=x[tet[t][2]];
        Vec3f s=x[tet[t][3]];
        std::cout<<t<<": "<<triple(q-p, r-p, s-p)<<std::endl;
    }
    */

    Gluvi::run();
    return 0;
}

