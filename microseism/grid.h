#ifndef GRID_H
#define GRID_H

/**
 * \file   grid.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Расчетная сетка.
 */

#include <vector>
#include <array>
#include <array>
#include <boost/multi_array.hpp>
#include "config.h"

/// Расчетная сетка.
class Grid {
    public:
	/// Конструктор.
	/**
	 * \param geometry   Имя файла с геометрией сетки.
	 * \param properties Имя файла со свойствами сетки.
	 */
	Grid(const std::string &geometry, const std::string &properties);

	Grid(const std::string &h5name);

	uint nx;		///< Число элементов вдоль оси X.
	uint ny;		///< Число элементов вдоль оси Y.
	uint nz;		///< Число элементов вдоль оси Z.

	std::vector<double> x;	///< Координаты узлов сетки по оси X.
	std::vector<double> y;	///< Координаты узлов сетки по оси Y.
	std::vector<double> z;	///< Координаты узлов сетки по оси Z.

	/// Свойства среды.
	struct Properties {
	    double young;	///< Модуль Юнга.
	    double poisson;	///< Коэффициент Пуассона.
	    double density;     ///< Плотность.
	    double stiffness;	///< Коэффициент при жесткости (<< 1).
	    double mass;	///< Коэффициент при массе (<< 1).
	    double viscosity;   ///< Коэффициент вязкости (1 - твердое тело, 0 - жидкость).
	};

	BoundingBox absorb;

	/// Свойства среды в ячейке с заданными координатами.
	/**
	 * \note Если в списке областей есть несколько подходящих, выдается
	 *       первая из них.
	 * \note Если в списке областей нет подходящей, выдается первая из
	 *       списка.
	 */
	Properties properties(uint i, uint j, uint k) const;

	/// Глобальный номер заданного узла заданного элемента.
	std::array<uint,8> node(uint i, uint j, uint k) const;

	/// Число узлов.
	uint num_nodes() const {
	    return nnodes;
	}

	/// Число элементов.
	inline uint num_elements() const {
	    return nelems;
	}

	/// Номер узла в глобальной нумерации.
	inline uint idx(uint i, uint j, uint k) const {
	    return i + j * nx_1 + k * nx_ny_1;
	}

	/// Размер элемента по оси x.
	inline real dx(uint i) const {
	    return x[i + 1] - x[i];
	}

	/// Размер элемента по оси y.
	inline real dy(uint i) const {
	    return y[i + 1] - y[i];
	}

	/// Размер элемента по оси z.
	inline real dz(uint i) const {
	    return z[i + 1] - z[i];
	}

	/// Принадлежность точки границе.
	inline bool boundary(std::array<uint,3> coord) const {
	    return coord[0] == 0 || coord[0] == nx ||
		   coord[1] == 0 || coord[1] == ny ||
		   coord[2] == nz;
	}

	/// Принадлежность точки фиксированному участку границы.
	inline bool boundary(uint i) const {
	    std::array<uint,3> p = {{
		i % (nx + 1),
		(i / (nx + 1)) % (ny + 1),
		i / ((nx + 1) * (ny + 1))
            }};

	    return boundary(p);
	}

    private:
	/// Область с заданными свойствами.
	struct Domain {
	    BoundingBox bbox;
	    Properties  prop;     ///< Свойства среды.
	};

	/// Области сетки с заданными свойствами.
	std::vector<Domain> domain;

        boost::multi_array<Properties, 3> P;

	uint nx_1;    // nx + 1
	uint nx_ny_1; // (nx + 1) * (ny + 1)
	uint nnodes;
	uint nelems;

        void fill_properties();
};

#endif
