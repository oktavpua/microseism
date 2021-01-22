#ifndef SAVER_H
#define SAVER_H

/**
 * \file   saver.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Сохранение сетки в файл.
 */

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <H5Cpp.h>
#include <pugixml.hpp>
#include "grid.h"
#include "partition.h"
#include "config.h"
#include "buildvoigt.h"

/// Сохранение сетки в формате XDMF/HDF5.
class HDFSaver {
    public:
	/// Конструктор.
	/**
	 * \param grid     Расчетная сетка.
	 * \param part     Разбиение сетки на процессы.
	 * \param mpi_rank Номер процесса.
	 * \param mpi_size Общее число процессов.
	 * \param fmask    Маска файлов для сохранения.
	 *
	 * \note Всего создается mpi_size+1 файлов. Мастер создает файл
	 * fmask.xmf, содержащий метаданные. Каждый процесс создает файл
	 * fmask-mpi_rank.h5, содержащий локальные данные процесса.
	 */
	HDFSaver(const Grid &grid, const Partition &part,
                const std::vector<int> &sp_ids,
                const std::vector<std::array<int,3>> &sp_coords,
		const std::string &fmask
		);

	/// Деструктор.
	~HDFSaver();

        /// Начальное условие.
        real read_checkpoint(std::vector<real> &q0, std::vector<real> &q1);

	/// Добавление очередного шага.
	void add_step(real time,
                const std::vector<real> &q0,
                const std::vector<real> &q1,
                const std::vector<real> &d
                );

	inline std::string file_name() const {
	    return hdfname[mpi_rank];
	}
    private:
	const Partition &part;
	std::vector<std::string> hdfname;
        std::vector<char>        do_savebox;
        std::vector<char>        do_savepnt;
        std::string xmf_name;

	H5::H5File hdf;
	H5::H5File hdf_chkpt;

	H5::FloatType         datatype;
	H5::DSetCreatPropList prop;
	H5::DataSpace         fbspace;
	H5::DataSpace         mbspace;
        H5::DataSpace         mpspace;
        H5::DataSpace         fpspace;
	H5::DataSet           time_ds;

        H5::DataSet           chkpt_time;
        H5::DataSet           chkpt_step;
        H5::DataSet           chkpt_q0;
        H5::DataSet           chkpt_q1;

        pugi::xml_document    xmf;
        pugi::xml_node        time_collection;

	std::vector<hsize_t> nx, ny, nz, np;
	hsize_t  step;

	BoundingBox sbox;

        std::chrono::high_resolution_clock::time_point last_checkpoint;
        bool have_chkpt;

        bool time_to_checkpoint() const;
        void save_checkpoint(
                real time,
                const std::vector<real> &q1,
                const std::vector<real> &d
                );
};

void save_matrix(
        const VoigtModel::Matrices &local,
        const VoigtModel::Matrices &remote,
        const std::vector<int>     &recv_nbr,
        const std::vector<uint>    &recv_idx,
        const std::vector<int>     &send_nbr,
        const std::vector<uint>    &send_idx,
        const std::vector<uint>    &send_col
        );

void read_matrix(
        VoigtModel::Matrices &local,
        VoigtModel::Matrices &remote,
        std::vector<int>     &recv_nbr,
        std::vector<uint>    &recv_idx,
        std::vector<int>     &send_nbr,
        std::vector<uint>    &send_idx,
        std::vector<uint>    &send_col
        );

#endif
