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
};

}

#endif
