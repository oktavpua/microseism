#ifndef BUILDVOIGT_H
#define BUILDVOIGT_H

/**
 * \file   buildvoigt.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Формирование матриц для модели Фойгта в упакованном виде.
 */

#include <vector>
#include <array>
#include <deque>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#include "grid.h"
#include "partition.h"
#include "config.h"

/// Формирование матриц для модели Фойгта в упакованном виде.
class VoigtModel {
    public:
	/// Набор матриц для модели Фойгта.
	/**
	 * Содержит матрицы, необходимые для продвижения по времени, в
	 * упакованном формате.
	 *
	 * Уникальные строки матриц хранятся в формате CSR. Вектор idx для
	 * каждой строки неупакованной матрицы содержит номер соответствующей
	 * ей строки в упакованной матрице.
	 */
	struct Matrices {
	    std::vector<uint>    idx;
	    std::vector<uint>    row;
	    std::vector<int> col;

	    std::vector<real>   M;
	    std::vector<real>   A1;
	    std::vector<real>   A2;

            size_t memory() const {
                return bytes(idx) +
                       bytes(row) +
                       bytes(col) +
                       bytes(M)   +
                       bytes(A1)  +
                       bytes(A2);
            }

            private:
                template <class Vector>
                static size_t bytes(const Vector &v) {
                    return v.size() * sizeof(v[0]);
                }
	};

	/// Конструктор.
	VoigtModel(const Grid &grid, const Partition &part,
		int mpi_rank, real tau,
		int first_plane, bool skip_first);

	/// Добавление значения в диагональную матрицу масс.
	void add(uint i, real m);

	/// Добавление значений в матрицы жесткости и демпфирования.
	void add(uint i, int j, real a1, real a2);

	/// Вызывается при переходе от одной плоскости x-y к другой.
	void shift_plane();

	/// Возвращает набор построенных упакованных матриц.
	Matrices matrices() const;

	/// Возвращает набор построенных упакованных матриц.
	Matrices remote_matrices() const;

	/// Число уникальных строк в матрицах.
	inline uint size() const { return unique_rows.size(); }
    private:
	const static uint row_width = 81;

	struct FullRow {
	    real M;

	    int len;

	    int col[row_width];
	    real      A1  [row_width];
	    real      A2  [row_width];

	    FullRow() : M(0), len(0) {}

	    inline void add(real m) { M += m; }

	    inline void add(int j, real a1, real a2) {
		if (len == 0) {
		    col[0] = j;
		    A1[0]  = a1;
		    A2[0]  = a2;

		    len++;
		} else {
		    int p = std::lower_bound(col, col + len, j) - col;

		    if (p < len && col[p] == j) {
			A1[p] += a1;
			A2[p] += a2;
		    } else {
			for(int q = len; q > p; q--) {
			    col[q] = col[q - 1];
			    A1[q]  = A1[q - 1];
			    A2[q]  = A2[q - 1];
			}

			col[p] = j;
			A1[p]  = a1;
			A2[p]  = a2;

			len++;
		    }
		}
	    }

	    inline void preprocess(real tau, bool bnd) {
		/* Уравнение движения волны имеет вид
		 *   M q_tt + A2 q_t + A1 q = f, (1)
		 * где M -- матрица масс, A1 -- матрица жесткости, f -- внешние силы, а
		 * матрица демпфирования A2 вычисляется как
		 *   A2 = alpha * A1 + beta * M.
		 * Уравнение (1) можно записать в виде
		 *   q_{i+1} = A_1 q_{i} + A_2 q_{i-1} + tau^2 * M^{-1} * f,
		 * где
		 *   A_1 = 2 - tau * (A2  + tau * A1) / M,
		 *   A_2 = tau * A2 / M - 1.
		 */
		real m = 1 / M;

		for(int i = 0; i < len; i++) {
		    A1[i] *= m;
		    A2[i] *= m;

		    if (bnd && col[i]) {
			A1[i] *= SCALE_BOUNDARY;
			A2[i] *= SCALE_BOUNDARY;
		    }
		}

		M = (bnd ? SCALE_BOUNDARY : 1) * tau * tau * m;
	    }
	};

	struct HashRow {
#ifdef DOUBLE_PRECISION
	    typedef uint mask_type;
#else
	    typedef unsigned int mask_type;
#endif
	    const static mask_type mask;
	    const static real      eps;

	    static std::hash<uint> hi;
	    static std::hash<real>   hd;

	    inline uint operator()(const FullRow *row) const {
		uint h = hi(row->len) ^ hd(trunc(row->M));

		for(int i = 0; i < row->len; i++) {
		    h ^= hi(row->col[i]);
		    h ^= hd(trunc(row->A1[i]));
		    h ^= hd(trunc(row->A2[i]));
		}

		return h;
	    }

	    inline real trunc(real v) const {
		if (fabs(v) < eps) return 0;

		union {
		    real      d;
		    mask_type i;
		} tmp;

		tmp.d  = v;
		tmp.i &= mask;

		return tmp.d;
	    }
	};

	struct CompRow {
	    const static real eps;

	    inline bool operator()(const FullRow *r1, const FullRow *r2) const {
		if (r1->len != r2->len)   return false;
		if (!equal(r1->M, r2->M)) return false;

		for(int i = 0; i < r1->len; i++) {
		    if (r1->col[i] != r2->col[i])   return false;
		    if (!equal(r1->A1[i], r2->A1[i])) return false;
		    if (!equal(r1->A2[i], r2->A2[i])) return false;
		}

		return true;
	    }

	    inline bool equal(real a, real b) const {
		return fabs(a - b) < eps;
	    }
	};

	const Grid      &grid;
	const Partition &part;

	int       mpi_rank;
	real      tau;
	uint plane_start, plane_size;
	bool      top_ready;

	std::vector<FullRow> top;
	std::vector<FullRow> bot;
	std::deque<FullRow>  unique_rows;

	std::vector<const FullRow*> ptr;

	std::unordered_set<const FullRow*, HashRow, CompRow> unique_idx;
};

#endif
