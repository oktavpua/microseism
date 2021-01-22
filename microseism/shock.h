#ifndef SHOCK_H
#define SHOCK_H

/**
 * \file   shock.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Информация об ударах.
 */

#include <vector>
#include <array>

#include "grid.h"
#include "partition.h"
#include "config.h"

/// Информация об ударах.
class Shock {
    public:
	/// Конструктор.
	/**
	 * \param mpi_rank Номер текущего процесса.
	 * \param grid     Расчетная сетка.
	 * \param part     Разбиение сетки на процессы.
	 * \param power    Имя файла с ударами.
	 */
	Shock(const Grid &grid, const Partition &part, const std::string &power);

	std::vector<double> time;	///< Моменты приложения силы.

	/// Целочисленные координаты приложения силы.
	std::vector<int> coord;

	/// Величины приложенных сил.
	std::vector<std::vector<std::array<double,3>>> F;
};

#endif
