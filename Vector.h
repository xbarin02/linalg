#include <vector>
#include <iostream>
#include <cmath>

namespace LinAlg {

template <typename T>
class Vector;

}

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Vector<T>& v);

namespace LinAlg {

template <typename T>
class Vector
{
private:
	// dimensions
	std::size_t size;
	std::vector<T> vector;
	Vector();
public:
	Vector(std::size_t size)
		: size(size)
	{
		vector.resize(size);
	}

	friend std::ostream& operator<< <>(std::ostream& os, const Vector& v);

	const T& operator[](std::size_t n) const
	{
		return vector[n];
	}

	T& operator[](std::size_t n)
	{
		return vector[n];
	}

	std::size_t dim() const
	{
		return size;
	}

	Vector operator+(const Vector &v)
	{
		*this += v;

		return *this;
	}

	Vector& operator+=(const Vector &v)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] += v.vector[n];
		}

		return *this;
	}

	Vector operator-(const Vector &v)
	{
		*this -= v;

		return *this;
	}

	Vector& operator-=(const Vector &v)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] -= v.vector[n];
		}

		return *this;
	}

	Vector operator*(T a)
	{
		*this *= a;

		return *this;
	}

	Vector& operator*=(T a)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] *= a;
		}

		return *this;
	}

	Vector operator/(T a)
	{
		*this /= a;

		return *this;
	}

	Vector& operator/=(T a)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] /= a;
		}

		return *this;
	}

	Vector operator-()
	{
		Vector v(*this);

		v *= -1;

		return v;
	}

	Vector operator+()
	{
		return *this;
	}

	// dot
	T operator*(const Vector &v) const
	{
		return *this *= v;
	}

	// dot
	T operator*=(const Vector &v) const
	{
		T s = 0;

		for (std::size_t n = 0; n < size; ++n) {
			s += vector[n] * v.vector[n];
		}

		return s;
	}

	// ||a|| = sqrt(a_i^2 + ... a_n^2) = sqrt(dot(a, a))
	// magnitude, L2-norm
	T length() const
	{
		T s = 0;

		for (std::size_t n = 0; n < size; ++n) {
			s += vector[n] * vector[n];
		}

		return std::sqrt(s);
	}
};

}

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Vector<T>& v)
{
	os << "Vector(";

	for (std::size_t n = 0; n < v.size; ++n) {
		os << v.vector[n];

		if (n != v.size - 1) {
			os << ", ";
		}
	}
	
	os << ")";

	return os;
}

template <typename T>
LinAlg::Vector<T> operator*(T a, const LinAlg::Vector<T>& v)
{
	LinAlg::Vector<T> w(v);

	return w * a;
}
