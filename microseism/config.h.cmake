#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <vector>
#include <string>

#include "util.h"
#include "bbox.h"
#include "profiler.h"

#cmakedefine DOUBLE_PRECISION
typedef @REAL@ real;
typedef unsigned int uint;

#define RTREE_BRANCHING  @RTREE_BRANCHING@
#define SCALE_BOUNDARY   @SCALE_BOUNDARY@
#define DOUBLE_HASH_MASK @DOUBLE_HASH_MASK@
#define DOUBLE_HASH_EPS  @DOUBLE_HASH_EPS@
#define FLOAT_HASH_MASK  @FLOAT_HASH_MASK@
#define FLOAT_HASH_EPS   @FLOAT_HASH_EPS@

#define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#ifdef DOUBLE_PRECISION
#  define MPI_MSREAL MPI_DOUBLE
extern unsigned long long hash_mask;
extern double hash_eps;
#else
#  define MPI_MSREAL MPI_FLOAT
extern unsigned int hash_mask;
extern float        hash_eps;
#endif

#define GITREVISION "${GITREVISION}"

extern int mpi_rank;
extern int mpi_size;

extern bool hdf_input;

extern std::string run_id;
extern std::vector<std::string> geometry_files;
extern std::string config_file;
extern std::string spoints_file;
extern std::string matrix_dir;
extern std::string chkpt_dir;
extern std::string fmask;
extern std::string done_file;
extern std::string scratch;

extern std::vector<std::string> force_files;
extern std::vector<std::string> force_ids;

extern real tmax;
extern real tau;
extern real wstep;
extern int  chkpt;
extern bool compress;

extern BoundingBox save_box;
extern BoundingBox absorb_box;

extern std::array<float,3> dim_weight;

extern profiler prof;

Log& log();

void read_config(int argc, char *argv[]);

// vim: ft=cpp
#endif
