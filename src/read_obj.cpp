#include "read_obj.h"
#include <fstream>
#include <cmath>
#include <cstdarg>
#include <cstdlib>

#define LINESIZE 1024 // maximum line size when reading .OBJ files

static bool
read_int(const char *s,
         int &value,
         bool &leading_slash,
         int &position)
{
    leading_slash=false;
    for(position=0; s[position]!=0; ++position){
        switch(s[position]){
            case '/':
                leading_slash=true;
                break;
            case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                goto found_int;
        }
    }
    return false;

found_int:
    value=0;
    for(;; ++position){
        switch(s[position]){
            case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
                value=10*value+s[position]-'0';
                break;
            default:
                return true;
        }
    }
    return true; // should never get here, but keeps compiler happy
}

static void
read_face_list(const char *s,
               std::vector<int> &vertex_list)
{
    vertex_list.clear();
    int v, skip;
    bool leading_slash;
    for(int i=0;;){
        if(read_int(s+i, v, leading_slash, skip)){
            if(!leading_slash)
                vertex_list.push_back(v-1); // correct for 1-based index
            //          vertex_list.push_back(v); // correct for 0-based index
            i+=skip;
        }else
            break;
    }
}

bool
read_objfile(std::vector<Vec3f>& x,
             std::vector<Vec3i>& tri,
             const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);
#ifdef WIN32
#define FILENAMELENGTH 256
   char filename[FILENAMELENGTH];
   _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
   std::ifstream input(filename, std::ifstream::binary);
#else
   char *filename;
   if (vasprintf(&filename, filename_format, ap) == -1) {
        return false;
   }
   std::ifstream input(filename, std::ifstream::binary);
   std::free(filename);
#endif
    va_end(ap);

    if(!input.good()) return false;

    x.clear();
    tri.clear();

    char line[LINESIZE];
    std::vector<int> vertex_list;
    while(input.good()){
        input.getline(line, LINESIZE);
        switch(line[0]){
            case 'v': // vertex data
                if(line[1]==' '){
                    Vec3f new_vertex;
                    std::sscanf(line+2, "%f %f %f", &new_vertex[0],
                                                    &new_vertex[1],
                                                    &new_vertex[2]);
                    x.push_back(new_vertex);
                }
                break;
            case 'f': // face data
                if(line[1]==' '){
                    read_face_list(line+2, vertex_list);
                    for(int j=0; j<(int)vertex_list.size()-2; ++j)
                        tri.push_back(Vec3i(vertex_list[0],
                                            vertex_list[j+1],
                                            vertex_list[j+2]));
                }
                break;
        }
    }
    return true;
}


bool
write_objfile(const std::vector<Vec3f>& x,
              const std::vector<Vec3i>& tri,
              const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);
#ifdef WIN32
#define FILENAMELENGTH 256
   char filename[FILENAMELENGTH];
   _vsnprintf(filename, FILENAMELENGTH, filename_format, ap);
   std::ofstream output(filename);
#else
   char *filename;
   if (vasprintf(&filename, filename_format, ap) == -1) {
       return false;
   }
   std::ofstream output(filename);
   std::free(filename);
#endif
    va_end(ap);
    
    if(!output.good()) return false;

    for (size_t v = 0; v < x.size(); ++v)
    {
        output << "v " << x[v][0] << " " << x[v][1] << " " << x[v][2] 
               << std::endl;
    }
    for (size_t f = 0; f < tri.size(); ++f)
    {
        output << "f " << tri[f][0]+1 << " " 
                       << tri[f][1]+1 << " " 
                       << tri[f][2]+1
               << std::endl;
    }
    
    return true;
}

