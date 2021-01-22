#include <iostream>
#include <iomanip>
#include <iterator>
#include <fstream>
#include <libconfig.h++>
#include <mpi.h>
#include <signal.h>
#include <vexcl/util.hpp>

#include <boost/scope_exit.hpp>

#include "config.h"
#include "util.h"
#include "grid.h"
#include "shock.h"
#include "partition.h"
#include "savepoints.h"
#include "msvoigt.h"

//---------------------------------------------------------------------------
void finalize(int) {
    if (mpi_rank == 0) {
	std::ofstream done(done_file);
	done << "Abnormal termination" << std::endl;
    }
    exit(100);
}

//---------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    // Обработка сигналов.
    signal(SIGINT,  finalize);
    signal(SIGTERM, finalize);

    // Инициализация MPI.
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    BOOST_SCOPE_EXIT(void) {
        MPI_Finalize();
    } BOOST_SCOPE_EXIT_END

    log() << "Version:    " << GITREVISION << Log::endl;
    log() << "World size: " << mpi_size << Log::endl;

    try {
        read_config(argc, argv);

	prof.tic("Read geometry");
	std::unique_ptr<Grid> grid;

        if (hdf_input) {
            log() << "Reading in \"" << geometry_files[0] << "\"" << Log::endl;
            grid.reset(new Grid(geometry_files[0]));
        } else {
            log() << "Reading in \"" << geometry_files[0] << "\" and \""
                  << geometry_files[1] << "\"" << Log::endl;
            grid.reset(new Grid(geometry_files[0], geometry_files[1]));
        }
        prof.toc("Read geometry");

        save_box.lo[0] = std::max<int>(save_box.lo[0], 0);
        save_box.lo[1] = std::max<int>(save_box.lo[1], 0);
        save_box.lo[2] = std::max<int>(save_box.lo[2], 0);
        save_box.hi[0] = std::min<int>(save_box.hi[0], grid->nx + 1);
        save_box.hi[1] = std::min<int>(save_box.hi[1], grid->ny + 1);
        save_box.hi[2] = std::min<int>(save_box.hi[2], grid->nz + 1);

        grid->absorb.lo[0] = std::max<int>(absorb_box.lo[0], 0);
        grid->absorb.lo[1] = std::max<int>(absorb_box.lo[1], 0);
        grid->absorb.lo[2] = std::max<int>(absorb_box.lo[2], 0);
        grid->absorb.hi[0] = std::min<int>(absorb_box.hi[0], grid->nx + 1);
        grid->absorb.hi[1] = std::min<int>(absorb_box.hi[1], grid->ny + 1);
        grid->absorb.hi[2] = std::min<int>(absorb_box.hi[2], grid->nz + 1);

	log() << "  tmax       = "  << tmax          << Log::endl
	      << "  tau        = "  << tau           << Log::endl
	      << "  wstep      = "  << wstep         << Log::endl
	      << "  chkpt      = "  << chkpt         << Log::endl
	      << "  fmask      = "  << fmask         << Log::endl
	      << "  scratch    = "  << scratch       << Log::endl
	      << "  compress   = "  << compress      << Log::endl
	      << "  hash_mask  = "  << std::hex      << hash_mask << std::dec << Log::endl
	      << "  hash_eps   = "  << hash_eps      << Log::endl
	      << "  savebox    = "  << save_box      << Log::endl
              << "  savepoints = "  << spoints_file  << Log::endl
	      << "  absorb     = "  << grid->absorb  << Log::endl
	      << "  dim_weight = [" << dim_weight[0] << ", "
	          		    << dim_weight[1] << ", "
	          		    << dim_weight[2] << "]"  << Log::endl;

	// Разбиение сетки на части.
	prof.tic("Grid partitioning");
	Partition part(*grid);
	log() << "Partitioning:\n" << part << Log::endl;
	prof.toc("Grid partitioning");

	// Чтение ударов.
	prof.tic("Read shocks");
        std::vector< Shock > shock;
        unsigned nshocks = force_files.size();
        shock.reserve(nshocks);
        for(auto &fname : force_files) {
            log() << "Reading in \"" << fname << "\"" << Log::endl;
            shock.emplace_back(*grid, part, fname);
        }
	prof.toc("Read shocks");

        // Чтение точек сохранения.
        std::vector< int > spoint_ids;
        std::vector< std::array<int, 3> > spoint_coords;
        if (!spoints_file.empty()) {
            read_save_points(spoints_file, spoint_ids, spoint_coords);
        }

	// Подготовка решателя.
	prof.tic("Solver setup");
	MSVoigt micros(*grid, part);
	prof.toc("Solver setup");

	// Решение.
	prof.tic("Solution");
	log() << "Advancing in time..." << Log::endl;
	micros.advance(spoint_ids, spoint_coords, shock);
	prof.toc("Solution");

	log() << prof;
    } catch (const std::exception &e) {
	log() << "Error: " << e.what() << Log::endl;

	if (mpi_rank == 0) {
	    std::ofstream done(done_file);
	    done << "Abnormal termination" << std::endl;
	    done << e.what() << std::endl;
	}

	return 1;
    }

    if (mpi_rank == 0) {
	std::ofstream done(done_file);
	done << "All done\n\n" << prof;
    }

    return 0;
}
