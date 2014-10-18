# Makefile for quartet
#
# This Makefile depends on a file called "Makefile.defs" which contains
# platform-specific definitions.
#
# To build quartet, first run "make depend" (which relies on g++ version 3.1 or
# higher) and then either "make" (equivalent to "make debug") or "make release".
#
# This is for GNU make; other versions of make may not run correctly.

MAIN_PROGRAM = quartet
SRC = predicates.cpp           \
      geometry_queries.cpp     \
      sdf.cpp                  \
      trimesh.cpp              \
      tet_mesh.cpp             \
      feature.cpp              \
      read_obj.cpp             \
      tet_quality.cpp          \
      match_features.cpp       \
      optimize_tet_mesh.cpp    \
      make_signed_distance.cpp \
      make_tet_mesh.cpp        \
      main.cpp

VIEWER_PROGRAM = view_tet
VIEWER_SRC = gluvi.cpp             \
             view_tet.cpp          \
             feature.cpp           \
             read_obj.cpp          \
             predicates.cpp        \
             geometry_queries.cpp  \
             tet_mesh.cpp          \
             tet_quality.cpp


SRC_ALL = $(SRC) $(VIEWER_SRC)

include Makefile.defs

# object files
RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC)))
VIEWER_RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(VIEWER_SRC)))
DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC)))
VIEWER_DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(VIEWER_SRC)))

.PHONY: all
all: debug

# how to make the main target (debug mode, the default)
$(MAIN_PROGRAM): $(DEBUG_OBJ)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS)

# how to make the main target (release mode)
$(MAIN_PROGRAM)_release: $(RELEASE_OBJ)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS)

# how to make the viewer application (debug mode, the default)
$(VIEWER_PROGRAM): $(VIEWER_DEBUG_OBJ)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS) $(VIEWER_LIBS)

# how to make the viewer application (release mode)
$(VIEWER_PROGRAM)_release: $(VIEWER_RELEASE_OBJ)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS) $(VIEWER_LIBS)

# how to compile the predicates.cpp source file
# This is different because optimization must be disabled 
# (at least, according to the TetGen source... not sure if it's true or not)
obj/predicates.o: src/predicates.cpp
	$(CC) -c $(RELEASE_FLAGS) -O0 -o $@ $^
obj_debug/predicates.o: src/predicates.cpp
	$(CC) -c $(DEBUG_FLAGS) -O0 -o $@ $^ 


.PHONY: release
release: $(MAIN_PROGRAM)_release $(VIEWER_PROGRAM)_release

.PHONY: debug
debug: $(MAIN_PROGRAM) $(VIEWER_PROGRAM)

# how to compile each file
.SUFFIXES:
obj/%.o:
	$(CC) -c $(RELEASE_FLAGS) -o $@ $<
obj_debug/%.o:
	$(CC) -c $(DEBUG_FLAGS) -o $@ $<

# cleaning up
.PHONY: clean
clean:
	-rm -f obj/*.o $(MAIN_PROGRAM) $(VIEWER_PROGRAM) obj_debug/*.o $(MAIN_PROGRAM)_release $(VIEWER_PROGRAM)_release *core

# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir obj
	-rm -f obj/depend
	$(foreach srcfile,$(SRC_ALL),$(DEPEND) -MM src/$(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) >> obj/depend;)
	-mkdir obj_debug
	-rm -f obj_debug/depend
	$(foreach srcfile,$(SRC_ALL),$(DEPEND) -MM src/$(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) >> obj_debug/depend;)

-include obj/depend
-include obj_debug/depend
