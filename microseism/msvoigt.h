#ifndef MSVOIGT_H
#define MSVOIGT_H

/**
 * \file   msvoigt.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Моделирование распространения сейсмических волн (модель Фойгта).
 */

#include <vector>
#include <array>
#include <set>
#include <string>
#include <mpi.h>

#include "util.h"
#include "grid.h"
#include "shock.h"
#include "partition.h"
#include "profiler.h"
#include "matrix.h"
#include "buildvoigt.h"

typedef unsigned int uint;

/// Моделирование распространения сейсмических волн.
class MSVoigt {
    public:
	/// Конструктор.
	MSVoigt(const Grid &grid, const Partition &part);

	/// Моделирование.
        void advance(
                const std::vector<int> &sp_ids,
                const std::vector<std::array<int,3>> &sp_coords,
		const std::vector<Shock> &shock
		);
    private:
	const Grid      &grid;		///< Расчетная сетка.
	const Partition &part;		///< Разбиение области.

	uint local_size;

        std::vector< std::vector<real> > q0;		///< Вектор перемещений.
	std::vector< std::vector<real> > q1;		///< Вектор перемещений.
	std::vector< std::vector<real> > q2;		///< Вектор перемещений.
	std::vector< std::vector<real> > q1_remote;	///< Вектор перемещений.
	std::vector< std::vector<real> > q2_remote;	///< Вектор перемещений.

	VoigtModel::Matrices mtx;
	VoigtModel::Matrices mtx_remote;

	struct {
	    std::vector<int>         nbr;
	    std::vector<uint>        idx;
	    std::vector<MPI_Request> req;
	} recv;

	struct {
	    std::vector<int>         nbr;
	    std::vector<uint>        idx;
	    std::vector<uint>        col;
            std::vector< std::vector<real> > val;
	    std::vector<MPI_Request> req;
	} send;

	void start_exchange(unsigned ns);
	void stop_exchange();
	void move_data(const std::string &scratch, const std::string &hdfname);
};

#endif
