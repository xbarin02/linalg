#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

namespace LinAlg {

template <typename T>
class Matrix
{
private:
	// columns
	std::size_t cols, rows;
	std::vector<T> column;
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

	void addColumnVector(const T &v)
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

	T getRowVector(std::size_t r) const
	{
		T v(cols);

		for (std::size_t c = 0; c < cols; ++c) {
			v[c] = column[c][r];
		}

		return v;
	}

	const T &getColumnVector(std::size_t c) const
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
			const T &v = getColumnVector(cc);
			// create subvector by eliminating r-th element
			T vv(rows - 1);
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

	typename T::type det() const
	{
		if (cols != rows) {
			throw std::domain_error("matrix must be square");
		}

		// end recursion
		if (cols == 1) {
			return column[0][0];
		}

		// expansion by 0-th row and c-th column
		typename T::type sum = 0;
		for (std::size_t c = 0; c < cols; ++c) {
			typename T::type a = getSubMatrix(0, c).det() * column[c][0];
			if (c % 2 == 0) {
				sum += a;
			} else {
				sum -= a;
			}
		}

		return sum;
	}
};

}

#endif
