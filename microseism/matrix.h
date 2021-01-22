#ifndef MATRIX_H
#define MATRIX_H

/**
 * \file   matrix.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Арифметические операции с небольшими плотными матрицами.
 */

#include <array>
#include "config.h"

/// Плотная матрица размера N x M.
template <int N, int M>
struct matrix {
    /// Число строк в матрице.
    static const int rows = N;

    /// Число столбцов в матрице.
    static const int cols = M;

    /// Матрица заполняется нулями.
    matrix() {
	for(int i = 0; i < N; i++)
	    for(int j = 0; j < M; j++)
		data[i][j] = 0;
    }

    /// Матрица заполняется значениями из списка инициализации.
#ifdef WIN32
    matrix(real s[N][M]) {
	for(int i = 0; i < N; i++)
	    for(int j = 0; j < M; j++)
		data[i][j] = s[i][j];
    }
#else
    matrix(std::initializer_list<real> s) {
	real *d = (real*)&data[0][0];

	for(auto v = s.begin(); v != s.end(); v++) *d++ = *v;
    }
#endif

    /// Добавление матрицы.
    matrix& operator+=(const matrix &A) {
	for(int i = 0; i < N; i++)
	    for(int j = 0; j < M; j++)
		data[i][j] += A(i,j);

	return *this;
    }

    /// Умножение на скаляр.
    matrix& operator*=(real a) {
	for(int i = 0; i < N; i++)
	    for(int j = 0; j < M; j++)
		data[i][j] *= a;

	return *this;
    }

    /// Элемент матрицы.
    real operator()(int i, int j) const {
	return data[i][j];
    }

    /// Элемент матрицы.
    real& operator()(int i, int j) {
	return data[i][j];
    }

    private:
	std::array<std::array<real,M>,N> data;
};

/// Транспонированная матрица.
/**
 * Хранит ссылку на исходную матрицу. Все операции сводятся к операциям с
 * исходной матрицей. Копирования данных не происходит.
 */
template <class Matrix>
struct TransposedMatrix {
    /// Число строк в матрице.
    static const int rows = Matrix::cols;

    /// Число столбцов в матрице.
    static const int cols = Matrix::rows;

    /// Конструктор.
    TransposedMatrix(const Matrix &M) : M(M) {}

    /// Элемент матрицы.
    real operator()(int i, int j) const {
	return M(j,i);
    }

    private:
	const Matrix &M;
};

/// Транспонирование матрицы.
template <class Matrix>
TransposedMatrix<Matrix> transp(const Matrix &A) {
    return TransposedMatrix<Matrix>(A);
}

/// Обращение матрицы 3x3.
inline matrix<3,3> inverse(const matrix<3,3> &A) {
#ifdef WIN32
    real buf[3][3] = {
	A(1,1) * A(2,2) - A(2,1) * A(1,2),
	A(0,2) * A(2,1) - A(2,2) * A(0,1),
	A(0,1) * A(1,2) - A(1,1) * A(0,2),
	A(1,2) * A(2,0) - A(2,2) * A(1,0),
	A(0,0) * A(2,2) - A(2,0) * A(0,2),
	A(0,2) * A(1,0) - A(1,2) * A(0,0),
	A(1,0) * A(2,1) - A(2,0) * A(1,1),
	A(0,1) * A(2,0) - A(2,1) * A(0,0),
	A(0,0) * A(1,1) - A(1,0) * A(0,1)
    };
    matrix<3,3> I = buf;
#else
    matrix<3,3> I = {
	A(1,1) * A(2,2) - A(2,1) * A(1,2),
	A(0,2) * A(2,1) - A(2,2) * A(0,1),
	A(0,1) * A(1,2) - A(1,1) * A(0,2),
	A(1,2) * A(2,0) - A(2,2) * A(1,0),
	A(0,0) * A(2,2) - A(2,0) * A(0,2),
	A(0,2) * A(1,0) - A(1,2) * A(0,0),
	A(1,0) * A(2,1) - A(2,0) * A(1,1),
	A(0,1) * A(2,0) - A(2,1) * A(0,0),
	A(0,0) * A(1,1) - A(1,0) * A(0,1)
    };
#endif

    real det = A(0,0) * I(0,0) + A(0,1) * I(1,0) + A(0,2) * I(2,0);

    I *= 1 / det;

    return I;
}

/// Умножение матриц.
template<class Left, class Right>
matrix<Left::rows, Right::cols> operator*(const Left &A, const Right &B) {
    static_assert(Left::cols == Right::rows, "Wrong matrix sizes");

    matrix<Left::rows, Right::cols> C;

    for(int i = 0; i < Left::rows; i++)
	for(int j = 0; j < Right::cols; j++)
	    for(int k = 0; k < Right::rows; k++)
		C(i,j) += A(i,k) * B(k,j);

    return C;
}

#endif
