#include <iostream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <mpi.h>

#include <boost/filesystem.hpp>

#include "msvoigt.h"
#include "saver.h"
#include "util.h"

#include "msvoigt_cl.h"

/*
 * [1] Голованов, Бережной. "Метод конечных элементов в механике деформируемых
 *     твердых тел"
 * [2] G. R. Liu, S. S. Quek. "The Finite Element Method. A Practical Course".
 */

enum MPI_tags {
    tagExcCols  = 100,
    tagExcVals  = 200,
    tagCombine0 = 300,
    tagCombine1 = 301,
};

//---------------------------------------------------------------------------
struct foreign_node {
    int    owner;
    uint id;

    foreign_node(int owner, uint id) : owner(owner), id(id) {}

    bool operator<(foreign_node node) const {
	return std::tie(owner, id) < std::tie(node.owner, node.id);
    }

    bool operator==(foreign_node node) const {
	return id == node.id;
    }
};

//---------------------------------------------------------------------------
MSVoigt::MSVoigt(
	const Grid &grid, const Partition &part
	)
    : grid(grid), part(part), local_size(3 * part.size(mpi_rank))
{
    if (boost::filesystem::exists(matrix_dir)) {
        prof.tic("Read setup");
        log() << "Reading setup from \"" << matrix_dir << "/\"..." << Log::endl;
        read_matrix(
                mtx, mtx_remote,
                recv.nbr, recv.idx,
                send.nbr, send.idx, send.col
                );
        if (mpi_size > 1) {
            recv.req.resize(recv.nbr.size());
            send.req.resize(send.nbr.size());
            send.val.resize(send.col.size());
        }
        prof.toc("Read setup");
    } else {
	log() << "Assembling matrices..." << Log::endl;
        BoundingBox bbox = part.box(mpi_rank);

        // Строим индекс:
        //   соответствие глобальных номеров узлов нашим локальным
        prof.tic("Indexing");
        std::unordered_map<uint,uint> g2l(2 * local_size);

        for(int k = bbox.lo[2], lid = 0; k < bbox.hi[2]; k++) {
            for(int j = bbox.lo[1]; j < bbox.hi[1]; j++) {
                for(int i = bbox.lo[0]; i < bbox.hi[0]; i++, lid++) {
                    uint gid = grid.idx(i, j, k);

                    for(int m = 0; m < 3; m++)
                        g2l[3 * gid + m] = 3 * lid + m;
                }
            }
        }
        prof.toc("Indexing");

        prof.tic("Assemble matrices");
        // Считаем dN / d_xi, dN / d_eta, dN / d_zeta в точках интегрирования.
        real corner[8][3] = {
            {-1, -1, -1},
            { 1, -1, -1},
            { 1,  1, -1},
            {-1,  1, -1},
            {-1, -1,  1},
            { 1, -1,  1},
            { 1,  1,  1},
            {-1,  1,  1}
        };

        // Производные функций формы по локальным координатам в узлах Гаусса.
        real gauss[] = {-1/sqrt(3.0), 1/sqrt(3.0)};

        real coord[3];
        matrix<3,8> dN[2][2][2];

        for(int k = 0; k < 2; k++) {
            coord[2] = gauss[k];

            for(int j = 0; j < 2; j++) {
                coord[1] = gauss[j];

                for(int i = 0; i < 2; i++) {
                    coord[0] = gauss[i];

                    for(int m = 0; m < 3; m++)
                        for(int n = 0; n < 8; n++)
                            dN[i][j][k](m,n) = 0.125 * corner[n][m] *
                                (1 + coord[(m + 1) % 3] * corner[n][(m + 1) % 3]) *
                                (1 + coord[(m + 2) % 3] * corner[n][(m + 2) % 3]);
                }
            }
        }

        // Формирование матриц масс, жесткости и демпфирования.
        /* Каждый процесс формирует свою часть глобальных матриц.
         * В множестве индексов remote_cols сохраняем номера интересующих нас
         * узлов, принадлежащих другим процессам.
         */
        std::set<foreign_node> remote_cols;

        uint elem_lo[] = {
            static_cast<uint>(bbox.lo[0] > 0 ? bbox.lo[0] - 1 : bbox.lo[0]),
            static_cast<uint>(bbox.lo[1] > 0 ? bbox.lo[1] - 1 : bbox.lo[1]),
            static_cast<uint>(bbox.lo[2] > 0 ? bbox.lo[2] - 1 : bbox.lo[2])
        };

        uint elem_hi[] = {
            std::min<uint>(bbox.hi[0], grid.nx),
            std::min<uint>(bbox.hi[1], grid.ny),
            std::min<uint>(bbox.hi[2], grid.nz)
        };

        std::unique_ptr<VoigtModel> lbuild(
                new VoigtModel(grid, part, mpi_rank, tau,
                    (int)elem_lo[2] - (int)bbox.lo[2], bbox.lo[2] > 0)
                );

        std::unique_ptr<VoigtModel> rbuild(
                new VoigtModel(grid, part, mpi_rank, tau,
                    (int)elem_lo[2] - (int)bbox.lo[2], bbox.lo[2] > 0)
                );

        for(uint k = elem_lo[2]; k < elem_hi[2]; k++) {
            for(uint j = elem_lo[1]; j < elem_hi[1]; j++) {
                for(uint i = elem_lo[0]; i < elem_hi[0]; i++) {
                    // Глобальные номера узлов элемента.
                    std::array<uint,8> node = grid.node(i, j, k);

                    // Свойства элемента.
                    auto prop = grid.properties(i, j, k);

                    std::array<uint,24> gid;
                    std::array<int,24>    owner;

                    for(int n = 0; n < 8; n++) {
                        gid[3 * n + 0] = 3 * node[n] + 0;
                        gid[3 * n + 1] = 3 * node[n] + 1;
                        gid[3 * n + 2] = 3 * node[n] + 2;

                        owner[3 * n + 0] = part[node[n]];
                        owner[3 * n + 1] = part[node[n]];
                        owner[3 * n + 2] = part[node[n]];
                    }

                    // Координаты углов элемента.
#ifdef WIN32
                    real corner_buf[8][3] = {
                        grid.x[i  ], grid.y[j  ], grid.z[k  ],
                        grid.x[i+1], grid.y[j  ], grid.z[k  ],
                        grid.x[i+1], grid.y[j+1], grid.z[k  ],
                        grid.x[i  ], grid.y[j+1], grid.z[k  ],
                        grid.x[i  ], grid.y[j  ], grid.z[k+1],
                        grid.x[i+1], grid.y[j  ], grid.z[k+1],
                        grid.x[i+1], grid.y[j+1], grid.z[k+1],
                        grid.x[i  ], grid.y[j+1], grid.z[k+1]
                    };
                    matrix<8,3> corner = corner_buf;
#else
                    matrix<8,3> corner = {
                        grid.x[i  ], grid.y[j  ], grid.z[k  ],
                        grid.x[i+1], grid.y[j  ], grid.z[k  ],
                        grid.x[i+1], grid.y[j+1], grid.z[k  ],
                        grid.x[i  ], grid.y[j+1], grid.z[k  ],
                        grid.x[i  ], grid.y[j  ], grid.z[k+1],
                        grid.x[i+1], grid.y[j  ], grid.z[k+1],
                        grid.x[i+1], grid.y[j+1], grid.z[k+1],
                        grid.x[i  ], grid.y[j+1], grid.z[k+1]
                    };
#endif

                    // Определитель Якобиана.
                    real detJ = 0.125 * grid.dx(i) * grid.dy(j) * grid.dz(k);

                    // На линейных элементах коэффициенты диагональной матрицы масс
                    // определяются путем деления массы всего элемента на число
                    // узлов элемента ([1], стр. 148).
                    real M_ii = prop.density * detJ;

                    {
                        std::array<real,2> v = {{
                            (2 - prop.mass * tau) * M_ii,
                            (prop.mass * tau - 1) * M_ii
                        }};

                        for(short n = 0; n < 8; n++) {
                            if (owner[n * 3] != mpi_rank) continue;

                            for(int m = 0; m < 3; m++) {
                                uint row = g2l[3 * node[n] + m];

                                lbuild->add(row, M_ii);
                                lbuild->add(row, 0, v[0], v[1]);
                                rbuild->add(row, M_ii);
                            }
                        }
                    }

                    // Матрица упругих постоянных.
                    real mu = prop.poisson;
                    real E  = prop.young;

#ifdef WIN32
                    real D_buf[6][6] = {
                        {1 - mu,   mu,     mu,       0,        0,        0},
                        {mu,   1 - mu,   mu,       0,        0,        0},
                        {mu,     mu,   1 - mu,     0,        0,        0},
                        {0,      0,      0,   0.5f - mu,    0,        0},
                        {0,      0,      0,       0,    0.5f - mu,    0},
                        {0,      0,      0,       0,        0,    0.5f - mu}
                    };
                    matrix<6,6> D = D_buf;
#else
                    matrix<6,6> D = {
                        1 - mu,   mu,     mu,       0,        0,        0,
                        mu,   1 - mu,   mu,       0,        0,        0,
                        mu,     mu,   1 - mu,     0,        0,        0,
                        0,      0,      0,   0.5f - mu,    0,        0,
                        0,      0,      0,       0,    0.5f - mu,    0,
                        0,      0,      0,       0,        0,    0.5f - mu
                    };
#endif

                    D *= E / ((1 + mu) * (1 - 2 * mu));

                    // Построение локальной матрицы жесткости.
                    /* Матрица жесткости элемента строится численным
                     * интегрированием следующей формулы:
                     *
                     *    K_e = \int_V [B^T] [D] [B] det(J) dV.
                     *
                     * Используется квадратура Гаусса.
                     */
                    matrix<24, 24> Ke;
                    for(int p = 0; p < 2; p++) {
                        for(int q = 0; q < 2; q++) {
                            for(int r = 0; r < 2; r++) {
                                matrix<3,8> dNxyz =
                                    inverse(dN[p][q][r] * corner) * dN[p][q][r];

                                matrix<6,24> B;

                                for(int n = 0; n < 8; n++) {
                                    B(0, 3 * n + 0) = dNxyz(0, n);
                                    B(1, 3 * n + 1) = dNxyz(1, n);
                                    B(2, 3 * n + 2) = dNxyz(2, n);
                                    B(3, 3 * n + 1) = dNxyz(2, n);
                                    B(3, 3 * n + 2) = dNxyz(1, n);
                                    B(4, 3 * n + 0) = dNxyz(2, n);
                                    B(4, 3 * n + 2) = dNxyz(0, n);
                                    B(5, 3 * n + 0) = dNxyz(1, n);
                                    B(5, 3 * n + 1) = dNxyz(0, n);
                                }

                                Ke += transp(B) * D * B;
                            }
                        }
                    }

                    Ke *= detJ;

                    real c1 = -tau * (prop.stiffness + tau);
                    real c2 = tau * prop.stiffness;

                    for(int n = 0; n < 24; n++) {
                        if (owner[n] != mpi_rank) continue;

                        uint row = g2l[gid[n]];

                        for(int m = 0; m < 24; m++) {
                            real k = Ke(n, m);

                            if (owner[m] == mpi_rank) {
                                lbuild->add(row, g2l[gid[m]] - row, c1 * k, c2 * k);
                            } else {
                                remote_cols.insert(foreign_node(owner[m], gid[m]));
                                rbuild->add(row, gid[m] - gid[n], c1 * k, c2 * k);
                            }
                        }
                    }
                }
            }
            lbuild->shift_plane();
            rbuild->shift_plane();

            log() << "." << Log::flush;

            if ((k - elem_lo[2]) % 50 == 49)
                log() << " " << 100 * (k - elem_lo[2]) / (elem_hi[2] - elem_lo[2])
                    << "% (" << lbuild->size() << ")" << Log::endl;
        }

        if (elem_hi[2] == grid.nz) {
            lbuild->shift_plane();
            rbuild->shift_plane();
        }

        log() << Log::endl;

        mtx = lbuild->matrices();

        log() << "Compression ratio: " << std::setprecision(2)
            << 100.0 * lbuild->size() / local_size << "%"
            << " (" << lbuild->size() << " unique rows)" << Log::endl;

        lbuild.reset();

        mtx_remote = rbuild->remote_matrices();
        rbuild.reset();

        {
            size_t bytes = mtx.memory() + mtx_remote.memory();
            size_t total;

            MPI_Allreduce(&bytes, &total, 1, MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

            log() << "Total matrix size: " << total / (1024.0 * 1024.0 * 1024.0) << " GB"
                << Log::endl;
        }

        prof.toc("Assemble matrices");

        // Настраиваем обмен между процессами
        prof.tic("Setup communication");
        std::vector<uint> recv_size(mpi_size, 0);
        std::vector<uint> recv_cols;
        std::unordered_map<uint,uint> r2l(2 * remote_cols.size());

        recv_cols.reserve(remote_cols.size());

        for(auto c = remote_cols.begin(); c != remote_cols.end(); c++) {
            recv_size[c->owner]++;
            r2l[c->id] = recv_cols.size();
            recv_cols.push_back(c->id);
        }

        std::vector<uint> comm_matrix(mpi_size * mpi_size);
        MPI_Allgather(recv_size.data(), mpi_size, MPI_UINT32_T,
                comm_matrix.data(), mpi_size, MPI_UINT32_T, MPI_COMM_WORLD);

        int rnbr  = 0;
        int snbr  = 0;
        int ssize = 0;

        for(int i = 0; i < mpi_size; i++) {
            if (recv_size[i]) {
                rnbr++;
            }

            if (uint nrecv = comm_matrix[mpi_size * i + mpi_rank]) {
                snbr++;
                ssize += nrecv;
            }
        }

        recv.nbr.resize(rnbr);
        recv.idx.resize(rnbr + 1);
        recv.req.resize(rnbr);

        send.nbr.resize(snbr);
        send.idx.resize(snbr + 1);
        send.col.resize(ssize);
        send.req.resize(snbr);

        recv.idx[0] = 0;
        for(int i = 0, j = 0; i < mpi_size; i++) {
            if (recv_size[i]) {
                recv.nbr[j++] = i;
                recv.idx[j]   = recv.idx[j - 1] + recv_size[i]; 
            }
        }

        send.idx[0] = 0;
        for(int i = 0, j = 0; i < mpi_size; i++) {
            if (uint nrecv = comm_matrix[mpi_size * i + mpi_rank]) {
                send.nbr[j++] = i;
                send.idx[j]   = send.idx[j - 1] + nrecv;
            }
        }

        for(uint i = 0; i < send.nbr.size(); i++)
            MPI_Irecv(&send.col[send.idx[i]], send.idx[i+1] - send.idx[i],
                    MPI_UINT32_T, send.nbr[i], tagExcCols, MPI_COMM_WORLD,
                    &send.req[i]);

        for(uint i = 0; i < recv.nbr.size(); i++)
            MPI_Isend(&recv_cols[recv.idx[i]], recv.idx[i+1] - recv.idx[i],
                    MPI_UINT32_T, recv.nbr[i], tagExcCols, MPI_COMM_WORLD,
                    &recv.req[i]);

        MPI_Waitall(recv.req.size(), recv.req.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(send.req.size(), send.req.data(), MPI_STATUSES_IGNORE);

        for(auto c = send.col.begin(); c != send.col.end(); c++)
            *c = g2l[*c];

        for(uint i = 0; i < local_size; i++)
            for(uint j = mtx_remote.row[i]; j < mtx_remote.row[i+1]; j++)
                mtx_remote.col[j] = r2l[mtx_remote.col[j]];
        prof.toc("Setup communication");

        {
            size_t bytes = (3 * local_size + 2 * recv.idx.back() + send.col.size()) * 8;
            size_t total;

            MPI_Allreduce(&bytes, &total, 1, MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

            log() << "Memory per force file: " << total / (1024.0 * 1024.0 * 1024.0) << " GB"
                << Log::endl;
        }

        prof.tic("Save setup info");
        log() << "Saving setup info..." << Log::endl;
        save_matrix(
                mtx, mtx_remote,
                recv.nbr, recv.idx,
                send.nbr, send.idx, send.col
                );
        prof.toc("Save setup info");
    }
}

//---------------------------------------------------------------------------
void MSVoigt::advance(
        const std::vector<int> &sp_ids,
        const std::vector<std::array<int,3>> &sp_coords,
        const std::vector<Shock> &shock
	)
{
    size_t ns = shock.size();

    tmax += 0.5 * tau;

    q0.resize(ns);
    q1.resize(ns);
    q2.resize(ns);
    q1_remote.resize(ns);

    send.val.resize(ns);

    for(unsigned s = 0; s < ns; ++s) {
        q0[s].resize(local_size, 0);
        q1[s].resize(local_size, 0);
        q2[s].resize(local_size, 0);
        q1_remote[s].resize(recv.idx.back());

        send.val[s].resize(send.col.size());
    }

    send.req.resize(send.nbr.size() * ns);
    recv.req.resize(recv.nbr.size() * ns);

    std::vector< std::unique_ptr<HDFSaver> > saver(ns);
    real start_time = 0;
    for(unsigned s = 0; s < ns; ++s) {
        saver[s].reset(
	    new HDFSaver(grid, part, sp_ids, sp_coords,
                fmask + "-" + force_ids[s]
                )
            );
        start_time = saver[s]->read_checkpoint(q0[s], q1[s]);
    }

    real next_chkpt = start_time + wstep;
    std::vector<uint> shock_num (ns, 0);
    for(unsigned s = 0; s < ns; ++s) {
        while (shock_num[s] < shock[s].time.size() && start_time > shock[s].time[shock_num[s]])
            shock_num[s]++;
    }

    MSVoigtCL solver(mtx.M, mtx.idx, mtx.row, mtx.col, mtx.A1, mtx.A2,
	    mtx_remote.row, mtx_remote.col, mtx_remote.A1, mtx_remote.A2,
	    send.col, recv.idx.back(), q0, q1);

    std::cout << mpi_rank << ": " << solver.device() << std::endl;

    // Цикл по времени.
    /* Вектор узловых смещений q обновляется в соответствии с уравнением
     *   M q_tt + C q_t + K q = 0
     */
    int iter = 0;
    for(real time = start_time; time <= tmax; iter++) {
	// Начинаем асинхронный обмен значениями с соседними процессами.
	solver.get_local_data(send.val);
	start_exchange(ns);

	// Вклад от локальных значений.
	solver.step_local();

	// Ждем окончания обмена.
	stop_exchange();
	solver.set_remote_data(q1_remote);

	// Вклад от значений с соседних процессов.
	solver.step_remote();

	// Задаем смещения из файла нагрузки.
        for(unsigned s = 0; s < ns; ++s) {
            if (shock_num[s] < shock[s].time.size() && time >= shock[s].time[shock_num[s]]) {
                log() << Log::endl << "Shock " << s << "/" << shock_num[s] << " at t = "
                    << time << Log::endl;

                solver.get_q0(s, q0[s]);

                for(uint i = 0; i < shock[s].coord.size(); i++) {
                    uint c = part.g2l(shock[s].coord[i]);

                    for(int m = 0; m < 3; m++)
                        q0[s][3 * c + m] += shock[s].F[shock_num[s]][i][m] * mtx.M[mtx.idx[3 * c + m]];
                }

                solver.set_q0(s, q0[s]);

                while (shock_num[s] < shock[s].time.size() && time >= shock[s].time[shock_num[s]])
                    shock_num[s]++;
            }
        }

        time += tau;

	// Запись в файл.
	if (time >= next_chkpt) {
	    prof.tic("Save data");
            solver.get_deriv(q2, tau);
            for(unsigned s = 0; s < ns; ++s) {
                solver.get_q0(s, q0[s]);
                solver.get_q1(s, q1[s]);
                saver[s]->add_step(time, q0[s], q1[s], q2[s]);
            }

	    while (time >= next_chkpt) next_chkpt += wstep;
	    prof.toc("Save data");
	}

	log() << "." << Log::flush;

	if (iter % 50 == 49)
	    log() << " " << int(100 * time / tmax) << "% (t = " << time << ")"
		<< Log::endl;
    }

    log() << Log::endl;

    // Перенос файлов из временной папки в текущую.
    if (scratch != ".") {
	prof.tic("Move data");
        for(unsigned s = 0; s < ns; ++s) {
            std::string hdfname = saver[s]->file_name();
            saver[s].reset(0);
            move_data(scratch, hdfname);
        }
	prof.toc("Move data");
    }
}

//---------------------------------------------------------------------------
void MSVoigt::start_exchange(unsigned ns) {
    for(unsigned s = 0, p = 0; s < ns; ++s)
        for(unsigned i = 0; i < recv.nbr.size(); i++, p++)
            MPI_Irecv(&q1_remote[s][recv.idx[i]], recv.idx[i+1] - recv.idx[i],
                    MPI_MSREAL, recv.nbr[i], tagExcVals, MPI_COMM_WORLD,
                    &recv.req[p]);

    for(unsigned s = 0, p = 0; s < ns; ++s)
        for(unsigned i = 0; i < send.nbr.size(); i++, p++)
            MPI_Isend(&send.val[s][send.idx[i]], send.idx[i+1] - send.idx[i],
                    MPI_MSREAL, send.nbr[i], tagExcVals, MPI_COMM_WORLD,
                    &send.req[p]);
}

//---------------------------------------------------------------------------
void MSVoigt::stop_exchange() {
    MPI_Waitall(recv.req.size(), recv.req.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(send.req.size(), send.req.data(), MPI_STATUSES_IGNORE);
}

//---------------------------------------------------------------------------
void MSVoigt::move_data(const std::string &scratch, const std::string &hdfname)
{
    log() << "Moving data to permanent location..." << Log::endl;
    for(int p = 0; p < mpi_size; p++) {
	if (p == mpi_rank && hdfname.size()) {
	    std::ostringstream cmd;
	    cmd << "mv " << scratch << "/" << hdfname << " " << hdfname;
	    if (0 != system(cmd.str().c_str())) {
		std::cout << "    failure!" << std::endl;
	    }
	}

	MPI_Barrier(MPI_COMM_WORLD);

	log() << "." << Log::flush;
	if (p % 50 == 49) {
	    log() << " " << 50 * p / mpi_size << "%" << Log::endl;
	}
    }
}
