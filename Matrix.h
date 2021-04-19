#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

namespace LinAlg {

template <typename T>
class Matrix;

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Matrix<T>& m);

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

	friend std::ostream& operator<< <>(std::ostream& os, const Matrix& m);

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

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Matrix<T>& m)
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

}

#endif
