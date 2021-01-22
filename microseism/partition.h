#ifndef PARTITION_H
#define PARTITION_H

/**
 * \file   partition.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Разбиение расчетной области на части.
 */

#include <iostream>
#include <vector>
#include <array>
#include "grid.h"
#include "rtree.h"

/// Разбиение расчетной области на части.
class Partition {
    public:
	/// Констуктор.
	/**
	 * \param grid     Расчетная область.
	 * \param mpi_size Число процессов (кол-во частей разбиения).
	 * \param mpi_rank Номер текущего процесса.
	 * \param wgt      Вес координат при разбиении. Координаты с большим
	 *		   весом бьются чаще.
	 */
	Partition(const Grid &grid);

	/// Номер процесса, владеющего указанным узлом сетки.
	/**
	 * \param i Номер узла.
	 */
	int operator[](uint i) const;

	/// Число узлов сетки, принадлежащих заданному процессу.
	/**
	 * \param p Номер процесса.
	 */
	uint size(int p) const;

	/// Прямоугольник, содержащий все узлы заданного процесса.
	/**
	 * \param p Номер процесса.
	 */
	BoundingBox box(int p) const;

	/// Преобразование глобального номера узла к локальному.
	/**
	 * \param idx Глобальный номер узла.
	 * \note Преобразование выполняется для текущего процесса.
	 */
	uint g2l(uint idx) const;

	/// Преобразование локального номера узла к глобальному.
	/**
	 * \param idx Локальный номер узла.
	 * \note Преобразование выполняется для текущего процесса.
	 */
	uint l2g(uint idx) const;

	friend std::ostream& operator<<(std::ostream &os, const Partition &part);
    private:
	std::vector<BoundingBox> bbox;
	RTree rtree;
	uint nx, ny, nz;

	void split(BoundingBox bb, int nparts);
};

#endif
