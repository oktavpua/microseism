#ifndef BBOX_H
#define BBOX_H

#include <array>
#include <limits>

/// Ограничивающий прямоугольник.
struct BoundingBox {
    /// Размерность.
    static const int dim = 3;

    /// Конструктор.
    /**
     * Создает бесконечный прямоугольник.
     */
    BoundingBox() {
	for(int i = 0; i < dim; i++) {
	    lo[i] = -std::numeric_limits<int>::max();
	    hi[i] =  std::numeric_limits<int>::max();
	}
    }

    /// Конструктор.
    /**
     * Создает прямоугольник, ограниченный координатами l и h.
     * \param lo Минимальные координаты.
     * \param hi Максимальные координаты.
     *
     * \note Координаты прямоугольника задаются полуоткрытым интервалом. Т.е.
     * для любой точки \f$x\f$ внутри прямоугольника \f$lo \leq x < hi\f$.
     */
    BoundingBox(
	    std::array<int,dim> lo,
	    std::array<int,dim> hi
	    ) : lo(lo), hi(hi) {}

    /// Конструктор.
    /**
     * Создает прямоугольник, содержащий единственную точку p.
     * \param p Координаты точки.
     */
    BoundingBox(std::array<int,dim> p) : lo(p), hi(p) {
	for (int i = 0; i < dim; i++) hi[i]++;
    }

    /// Проверка на пересечение.
    bool  overlaps(const BoundingBox &bb) const {
	for(int i = 0; i < dim; i++)
	    if (lo[i] >= bb.hi[i] || bb.lo[i] >= hi[i]) return false;
	return true;
    }

    bool owns(std::array<int, dim> p) const {
        for(int i = 0; i < dim; ++i)
            if (lo[i] > p[i] || p[i] >= hi[i]) return false;
        return true;
    }

    /// Площадь пересечения.
    int overlap(const BoundingBox &bb) const {
	int area = 1;
	for (int i = 0; area && i < dim; i++) {
	    // this makes it easier to understand
	    const int x1 = lo[i];
	    const int x2 = hi[i];
	    const int y1 = bb.lo[i];
	    const int y2 = bb.hi[i];

	    // left edge outside left edge
	    if (x1 <= y1) {
		// and right edge inside left edge
		if (y1 < x2) {
		    // right edge outside right edge
		    if (y2 < x2)
			area *= y2 - y1;
		    else
			area *= x2 - y1;

		    continue;
		}
	    }

	    // right edge inside left edge
	    else if (x1 < y2) {
		// right edge outside right edge
		if (x2 < y2)
		    area *= x2 - x1;
		else
		    area *= y2 - x1;

		continue;
	    }

	    // if we get here, there is no overlap
	    return 0;
	}

	return area;
    }

    /// Объединение.
    void  stretch(const BoundingBox &bb) {
	for(int i = 0; i < dim; i++) {
	    if (lo[i] > bb.lo[i]) lo[i] = bb.lo[i];
	    if (hi[i] < bb.hi[i]) hi[i] = bb.hi[i];
	}
    }

    /// Площадь (объем).
    int area() const {
	int a = 1;
	for(int i = 0; i < dim; i++)
	    a *= hi[i] - lo[i];
	return a;
    }

    std::array<int,dim> lo; ///< Точка с минимальными координатами.
    std::array<int,dim> hi; ///< Тkочка с максимальными координатами.
};

std::ostream& operator<<(std::ostream &os, const BoundingBox &bb);



#endif
