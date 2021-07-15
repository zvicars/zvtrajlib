# Prefer to keep flags private)
add_library (XDRFILE
trr_seek.h
trr_seek.c
xdrfile_trr.h
xdrfile_trr.c
xdrfile_xtc.h
xdrfile_xtc.c
xdrfile.h
xdrfile.c
xtc_seek.h
xtc_seek.c
)
# Add source files
target_include_directories(XDRFILE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})