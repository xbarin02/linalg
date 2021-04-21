#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <utility>
#include "Vector.h"

namespace LinAlg {

template <typename T>
class Matrix
{
private:
	// columns
	std::size_t cols, rows;
	std::vector<Vector<T>> column;
public:
	Matrix()
		: cols(0), rows(0), column(0)
	{}

	Matrix(const Matrix &v)
		: cols(v.cols), rows(v.rows), column(v.column)
	{}

	friend void swap(Matrix &first, Matrix &second)
	{
		using std::swap;

		swap(first.cols, second.cols);
		swap(first.rows, second.rows);
		swap(first.column, second.column);
	}

	Matrix& operator=(Matrix v)
	{
		swap(*this, v);

		return *this;
	}

	friend std::ostream& operator<<(std::ostream& os, const Matrix& m)
	{
		for (std::size_t r = 0; r < m.rows; ++r) {
			os << "[ ";

			for (std::size_t c = 0; c < m.cols; ++c) {

				os << m.column[c][r];

				if (c + 1 != m.cols) {
					os << ", ";
				}
			}

			os << " ]\n";
		}

		return os;
	}

	void addColumnVector(const Vector<T> &v)
	{
		if (cols == 0) {
			// add anything
			column.push_back(v);
			cols = 1;
			rows = v.dim();
		} else {
			// check whether number of rows match
			if (rows == v.dim()) {
				column.push_back(v);
				cols++;
			} else {
				throw std::domain_error("wrong number of rows");
			}
		}
	}

	Vector<T> getRowVector(std::size_t r) const
	{
		Vector<T> v(cols);

		for (std::size_t c = 0; c < cols; ++c) {
			v[c] = column[c][r];
		}

		return v;
	}

	void swapRowVectors(std::size_t r0, std::size_t r1)
	{
		for (std::size_t c = 0; c < cols; ++c) {
			std::swap(column[c][r0], column[c][r1]);
		}
	}

	// add row vector 'v' to r-th row in the matrix
	void addVectorToRow(std::size_t r, const Vector<T> &v)
	{
		if (cols != v.dim()) {
			throw std::domain_error("wrong number of elements");
		}

		for (std::size_t c = 0; c < cols; ++c) {
			column[c][r] += v[c];
		}
	}

	void subVectorFromRow(std::size_t r, const Vector<T> &v)
	{
		if (cols != v.dim()) {
			throw std::domain_error("wrong number of elements");
		}

		for (std::size_t c = 0; c < cols; ++c) {
			column[c][r] -= v[c];
		}
	}

	void mulRowVector(std::size_t r, T factor)
	{
		for (std::size_t c = 0; c < cols; ++c) {
			column[c][r] *= factor;
		}
	}

	// https://en.wikipedia.org/wiki/Row_echelon_form#Reduced_row_echelon_form
	void toRref()
	{
		std::size_t lead = 0;

		for (std::size_t r = 0; r < rows; ++r) {
			if (cols <= lead) {
				return;
			}
			std::size_t i = r;
			while (column[lead][i] == 0) {
				i++;
				if (rows == i) {
					i = r;
					lead++;
					if (cols == lead) {
						return;
					}
				}
			} /* end while */
			if (i != r) {
				swapRowVectors(i, r);
			}
			T factor = column[lead][r];
			mulRowVector(r, 1/factor); // FIXME: explicitly set column[lead][r] = 11
			for (i = 0; i < rows; ++i) {
				if (i != r) {
					T factor = column[lead][i];
					subVectorFromRow(i, getRowVector(r) * factor); // FIXME: set zeros?
				}
			} /* end for */
			lead++;
		} /* end for */
	}

	const Vector<T> &getColumnVector(std::size_t c) const
	{
		return column[c];
	}

	Vector<T> &getColumnVector(std::size_t c)
	{
		return column[c];
	}

	Matrix getTranspose() const
	{
		Matrix m;

		for (std::size_t r = 0; r < rows; ++r) {
			m.addColumnVector(getRowVector(r));
		}

		return m;
	}

	// returns submatrix created by eliminating r-th row and c-th column
	Matrix getSubMatrix(std::size_t r, std::size_t c) const
	{
		Matrix m;

		// for each column vector except c-th one
		for (std::size_t cc = 0; cc < cols; ++cc) {
			if (cc == c) {
				continue;
			}
			// getColumnVector()
			const Vector<T> &v = getColumnVector(cc);
			// create subvector by eliminating r-th element
			Vector<T> vv(rows - 1);
			std::size_t i = 0;
			for (std::size_t rr = 0; rr < rows; ++rr) {
				if (rr == r) {
					continue;
				}
				vv[i++] = v[rr];
			}
			// addColumnVector()
			m.addColumnVector(vv);
		}

		return m;
	}

	T det() const
	{
		if (cols != rows) {
			throw std::domain_error("matrix must be square");
		}

		// end recursion
		if (cols == 1) {
			return column[0][0];
		}

		// expansion by 0-th row and c-th column
		T sum = 0;
		for (std::size_t c = 0; c < cols; ++c) {
			T a = getSubMatrix(0, c).det() * column[c][0];
			if (c % 2 == 0) {
				sum += a;
			} else {
				sum -= a;
			}
		}

		return sum;
	}

	// matrix-matrix product
	Matrix operator*(const Matrix &right) const
	{
		Matrix result;
		Matrix leftTransposed = getTranspose();

		for (std::size_t c = 0; c < right.cols; ++c) {
			Vector<T> v(cols);

			for (std::size_t r = 0; r < cols; ++r) {
				v[r] = leftTransposed.getColumnVector(r) * right.getColumnVector(c);
			}

			result.addColumnVector(v);
		}

		return result;
	}

	Matrix& operator*=(const Matrix &m)
	{
		*this = *this * m;

		return *this;
	}

	Matrix operator+(const Matrix &v) const
	{
		Matrix copy(*this);

		copy += v;

		return copy;
	}

	Matrix& operator+=(const Matrix &m)
	{
		if (cols != m.cols || rows != m.rows) {
			throw std::domain_error("matrix dimensions must agree");
		}

		for (std::size_t c = 0; c < cols; ++c) {
			getColumnVector(c) += m.getColumnVector(c);
		}

		return *this;
	}

	static Matrix identity(std::size_t dim)
	{
		Vector<T> v(dim);

		for (std::size_t n = 0; n < dim; ++n) {
			v[n] = 0;
		}

		Matrix m;

		for (std::size_t n = 0; n < dim; ++n) {
			v[n] = 1;
			m.addColumnVector(v);
			v[n] = 0;
		}

		return m;
	}
};

}

#endif
