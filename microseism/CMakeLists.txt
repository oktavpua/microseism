cmake_minimum_required(VERSION 2.6)
project(microseism)

if (WIN32)
    set(Boost_USE_STATIC_LIBS ON)
endif (WIN32)
find_package(Boost COMPONENTS
    chrono
    date_time
    filesystem
    program_options
    regex
    system
    )

include_directories(${Boost_INCLUDE_DIRS})

execute_process(
    COMMAND git rev-parse --short HEAD
    OUTPUT_VARIABLE GITREVISION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

set(REAL            double CACHE STRING "Type of real numbers (double/float)")
set(RTREE_BRANCHING 8      CACHE STRING "R-Tree maximum branching")

set(SCALE_BOUNDARY   1e-6               CACHE STRING "Inverted scale for boundary values")
set(DOUBLE_HASH_MASK 0xffffffffff000000 CACHE STRING "Hash mask for double")
set(FLOAT_HASH_MASK  0xfffffc00         CACHE STRING "Hash mask for float")
set(DOUBLE_HASH_EPS  1e-5               CACHE STRING "Hash eps for double")
set(FLOAT_HASH_EPS   1e-4               CACHE STRING "Hash eps for float")

if (REAL STREQUAL "double")
    set(DOUBLE_PRECISION ON)
endif (REAL STREQUAL "double")

if (WIN32)
    set (HDFROOT "C:/HDF5" CACHE STRING "HDF5 root folder")
    set (MPIROOT "C:/MPICH2" CACHE STRING "MPICH2 root folder")

    find_library(LIBCONFIG libconfig++ c:/lib)
    find_path(LIBCONFIG_INCLUDE libconfig.h++ c:/lib)

    find_library(LIBHDF5 hdf5 ${HDFROOT}/lib)
    find_library(LIBHDF5_CPP hdf5_cpp ${HDFROOT}/lib)
    find_path(LIBHDF5_INCLUDE H5Cpp.h ${HDFROOT}/include)

    find_library(MPICH2_LIBRARY mpi ${MPIROOT}/lib)
    find_path(MPICH2_INCLUDE mpi.h ${MPIROOT}/include)

    include_directories(
	${LIBCONFIG_INCLUDE}
	${LIBHDF5_INCLUDE}
	${MPICH2_INCLUDE}
	)
else (WIN32)
    find_library(LIBCONFIG config++ $ENV{HOME}/lib)
    find_path(LIBCONFIG_INCLUDE libconfig.h++ $ENV{HOME}/include)

    include_directories(${LIBCONFIG_INCLUDE})
endif (WIN32)

set(SRC
    buildvoigt.cpp
    config.cpp
    grid.cpp
    msvoigt_cl.cpp
    msvoigt.cpp
    partition.cpp
    profiler.cpp
    rtree.cpp
    savepoints.cpp
    saver.cpp
    shock.cpp
    util.cpp
    )

#----------------------------------------------------------------------------
# Find VexCL
#----------------------------------------------------------------------------
set(VEXCL_ROOT $ENV{VEXCL_ROOT} CACHE STRING "VexCL root")
if (VEXCL_ROOT)
    include_directories( ${VEXCL_ROOT} )

    set(VEXCL_BACKEND "OpenCL" CACHE STRING "Select VexCL backend (OpenCL/CUDA)")
    set_property(CACHE VEXCL_BACKEND PROPERTY STRINGS "OpenCL" "CUDA")

    #----------------------------------------------------------------------------
    # Find VexCL backend
    #----------------------------------------------------------------------------
    if ("${VEXCL_BACKEND}" STREQUAL "OpenCL")
        set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
        find_package(OpenCL REQUIRED)
        include_directories( ${OPENCL_INCLUDE_DIRS} )
        set(VEXCL_BACKEND_LIBRARIES ${OPENCL_LIBRARIES})
        add_definitions(-DVEXCL_BACKEND_OPENCL)
    elseif ("${VEXCL_BACKEND}" STREQUAL "CUDA")
        find_package(CUDA REQUIRED)
        include_directories( ${CUDA_INCLUDE_DIRS} )
        set(VEXCL_BACKEND_LIBRARIES ${CUDA_CUDA_LIBRARY})
        add_definitions(-DVEXCL_BACKEND_CUDA)
    endif()
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/config.h
    )

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${LIBCONFIG_INCLUDE}
    )

add_executable(microseism driver.cpp ${SRC})
add_executable(microseism_setup setup.cpp ${SRC})
add_executable(convertprop convertprop.cpp grid.cpp rtree.cpp util.cpp)
add_executable(txt2hdf txt2hdf.cpp util.cpp)

target_link_libraries(microseism ${VEXCL_BACKEND_LIBRARIES})
target_link_libraries(microseism_setup ${VEXCL_BACKEND_LIBRARIES})

if (WIN32)
    target_link_libraries(microseism
	${LIBCONFIG}
	${LIBHDF5}
	${LIBHDF5_CPP}
	${MPICH2_LIBRARY}
        ${Boost_LIBRARIES}
	)
    target_link_libraries(microseism_setup
	${LIBCONFIG}
	${LIBHDF5}
	${LIBHDF5_CPP}
	${MPICH2_LIBRARY}
        ${Boost_LIBRARIES}
	)
    target_link_libraries(convertprop
	${LIBCONFIG}
	${LIBHDF5}
	${LIBHDF5_CPP}
	)
    target_link_libraries(txt2hdf
	${LIBHDF5}
	${LIBHDF5_CPP}
	)
else (WIN32)
    target_link_libraries(microseism
	${LIBCONFIG}
        ${Boost_LIBRARIES}
	hdf5 hdf5_cpp pugixml
	)
    target_link_libraries(microseism_setup
	${LIBCONFIG}
        ${Boost_LIBRARIES}
	hdf5 hdf5_cpp pugixml
	)
    target_link_libraries(convertprop
	${LIBCONFIG}
	hdf5 hdf5_cpp
	)
    target_link_libraries(txt2hdf
	hdf5 hdf5_cpp
	)

    set_target_properties(microseism
	PROPERTIES
        COMPILE_FLAGS "-std=c++0x -Wall -Wno-comment -Wno-ignored-attributes -Wno-sign-compare"
	)
    set_target_properties(microseism_setup
	PROPERTIES
        COMPILE_FLAGS "-std=c++0x -Wall -Wno-comment -Wno-ignored-attributes -Wno-sign-compare"
	)
    set_target_properties(convertprop
	PROPERTIES
	COMPILE_FLAGS "-std=c++0x"
	)
    set_target_properties(txt2hdf
	PROPERTIES
	COMPILE_FLAGS "-std=c++0x"
	)

    add_custom_target(tags
	COMMAND ctags --exclude='build/*' -R ..
	)
endif (WIN32)

add_subdirectory(doc)
add_subdirectory(gendata)

install(TARGETS microseism       RUNTIME DESTINATION bin)
install(TARGETS microseism_setup RUNTIME DESTINATION bin)
