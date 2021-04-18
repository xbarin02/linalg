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
