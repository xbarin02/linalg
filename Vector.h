#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <iostream>
#include <cmath>

namespace LinAlg {

template <typename T>
class Vector;

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Vector<T>& v);

template <typename T>
class Vector
{
private:
	// dimensions
	std::size_t size;
	std::vector<T> vector;
public:
	Vector()
		: size(0), vector(0)
	{}

	Vector(std::size_t size)
		: size(size), vector(size)
	{}

	Vector(const Vector &v)
		: size(v.size), vector(v.vector)
	{}

	Vector(std::initializer_list<T> l)
		: size(l.size()), vector(l)
	{}

	friend void swap(Vector &first, Vector &second)
	{
		using std::swap;

		swap(first.size, second.size);
		swap(first.vector, second.vector);
	}

	Vector& operator=(Vector v)
	{
		swap(*this, v);

		return *this;
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

	Vector operator+(const Vector &v) const
	{
		Vector copy(*this);

		copy += v;

		return copy;
	}

	Vector& operator+=(const Vector &v)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] += v.vector[n];
		}

		return *this;
	}

	Vector operator-(const Vector &v) const
	{
		Vector copy(*this);

		copy -= v;

		return copy;
	}

	Vector& operator-=(const Vector &v)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] -= v.vector[n];
		}

		return *this;
	}

	Vector operator*(T a) const
	{
		Vector copy(*this);

		copy *= a;

		return copy;
	}

	Vector& operator*=(T a)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] *= a;
		}

		return *this;
	}

	Vector operator/(T a) const
	{
		Vector copy(*this);

		copy /= a;

		return copy;
	}

	Vector& operator/=(T a)
	{
		for (std::size_t n = 0; n < size; ++n) {
			vector[n] /= a;
		}

		return *this;
	}

	Vector operator-() const
	{
		Vector copy(*this);

		copy *= -1;

		return copy;
	}

	Vector operator+() const
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

	std::size_t dim() const
	{
		return size;
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

	T angle(const Vector &v) const
	{
		const Vector &a = *this;
		const Vector &b = v;

		return std::acos((a * b) / a.length() / b.length());
	}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const LinAlg::Vector<T>& v)
{
	os << "Vector(";

	for (std::size_t n = 0; n < v.size; ++n) {
		os << v.vector[n];

		if (n + 1 != v.size) {
			os << ", ";
		}
	}
	
	os << ")";

	return os;
}

}

template <typename T>
LinAlg::Vector<T> operator*(T a, const LinAlg::Vector<T>& v)
{
	LinAlg::Vector<T> w(v);

	return w * a;
}

#endif
