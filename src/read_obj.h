#ifndef READ_OBJ_H
#define READ_OBJ_H

#include <vector>
#include "vec.h"

bool
read_objfile(std::vector<Vec3f>& x,
             std::vector<Vec3i>& tri,
             const char *filename_format, ...);

bool
write_objfile(const std::vector<Vec3f>& x,
              const std::vector<Vec3i>& tri,
              const char *filename_format, ...);

#endif
