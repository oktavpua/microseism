#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <libconfig.h++>

#include "config.h"

int mpi_rank;
int mpi_size;

#ifdef DOUBLE_PRECISION
unsigned long long hash_mask = DOUBLE_HASH_MASK;
double             hash_eps  = DOUBLE_HASH_EPS;
#else
unsigned int hash_mask = FLOAT_HASH_MASK;
float        hash_eps  = FLOAT_HASH_EPS;
#endif

bool hdf_input;

std::string run_id;
std::vector<std::string> geometry_files;
std::string config_file;
std::string spoints_file;
std::string matrix_dir;
std::string chkpt_dir;
std::string done_file;
std::string fmask   = "microseism";
std::string scratch = ".";

std::vector<std::string> force_files;
std::vector<std::string> force_ids;

profiler prof;

real tmax       = 10;
real tau        = 1e-3;
real wstep      = 1e-1;
int  chkpt      = 10;
bool compress   = true;

BoundingBox save_box;
BoundingBox absorb_box;

std::array<float,3> dim_weight = {{1, 1, 1}};

Log& log() {
    static Log l(std::cout, mpi_rank);
    return l;
}

void read_config(int argc, char *argv[]) {
    namespace po = boost::program_options;
    namespace fs = boost::filesystem;
    po::options_description desc("Options");

    std::string uuid;

    desc.add_options()
        ("help,h", "Show help")
        (
         "hdf_input,b",
         po::bool_switch(&hdf_input),
         "Input files are in binary HDF5 format"
        )
        (
         "geometry,g",
         po::value<std::vector<std::string>>(&geometry_files)->required()->multitoken(),
         "Geometry and medium properties file name(s). "
         "In case of binary input (hdf_input is set), this is single hdf file."
         "Otherwise, two text files should be provided "
         "(geometry and medium properties)."
        )
        (
         "config,c",
         po::value<std::string>(&config_file)->required(),
         "Configuration file name"
        )
        (
         "force,f",
         po::value<std::vector<std::string>>(&force_files)->required()->multitoken(),
         "Force file names (one or more)"
        )
        (
         "uuid,i",
         po::value<std::string>(&uuid),
         "Unique job id. Used as a postfix for 'done' file and 'matrix' folder"
        )
        (
         "matrix,m",
         po::value<std::string>(&matrix_dir),
         "Directory with precomputed matrix for the problem"
        )
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        if (mpi_rank == 0)
            std::cout << desc << std::endl;
        exit(0);
    }

    po::notify(vm);

    if (hdf_input)
        precondition(
            geometry_files.size() == 1,
            "hdf_input is set: "
            "geometry should be specified as a single hdf file."
            );
    else
        precondition(
            geometry_files.size() == 2,
            "hdf_input is not set: "
            "geometry should be specified as two files "
            "(geometry and medium properties)."
            );

    if (!uuid.empty()) uuid = "-" + uuid;

    if (matrix_dir.empty()) matrix_dir = "matrix" + uuid;

    chkpt_dir = "chkpt" + uuid;
    done_file = "done" + uuid + ".txt";

    for(const auto &f : force_files)
        force_ids.push_back(fs::basename(f));

    log() << "Parsing \"" << config_file << "\"" << Log::endl;

    try {
        libconfig::Config cfg;
        cfg.setAutoConvert(true);
        cfg.readFile(config_file.c_str());

        cfg.lookupValue("tmax",       tmax);
        cfg.lookupValue("tau",        tau);
        cfg.lookupValue("wstep",      wstep);
        cfg.lookupValue("chkpt",      chkpt);
        cfg.lookupValue("fmask",      fmask);
        cfg.lookupValue("compress",   compress);
        cfg.lookupValue("scratch",    scratch);
        hash_mask = ~hash_mask;
        cfg.lookupValue("hash_mask",  hash_mask);
        hash_mask = ~hash_mask;
        cfg.lookupValue("hash_eps",   hash_eps);

        if(cfg.exists("savebox")) {
            libconfig::Setting &bb = cfg.lookup("savebox");

            save_box.lo[0] = bb[0];
            save_box.lo[1] = bb[1];
            save_box.lo[2] = bb[2];
            save_box.hi[0] = bb[3];
            save_box.hi[1] = bb[4];
            save_box.hi[2] = bb[5];
        }

        cfg.lookupValue("savepoints",   spoints_file);

        if(cfg.exists("absorb")) {
            libconfig::Setting &bb = cfg.lookup("absorb");

            absorb_box.lo[0] = bb[0];
            absorb_box.lo[1] = bb[1];
            absorb_box.lo[2] = bb[2];
            absorb_box.hi[0] = bb[3];
            absorb_box.hi[1] = bb[4];
            absorb_box.hi[2] = bb[5];
        }

        if (cfg.exists("dim_weight")) {
            libconfig::Setting &dw = cfg.lookup("dim_weight");

            dim_weight[0] = dw[0];
            dim_weight[1] = dw[1];
            dim_weight[2] = dw[2];
        }
    } catch(const libconfig::FileIOException&) {
        precondition(false, "Failed to read configuration file");
    }
}

