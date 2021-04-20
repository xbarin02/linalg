#include "Vector.h"
#include "Matrix.h"

#include <ostream>

using namespace LinAlg;

int main()
{
	Vector<float> v(5);

	std::cout << v << "\n";

	v[2] = 3.14;
	std::cout << v[2] << "\n";

	std::cout << v << "\n";

	Vector<float> w;
	std::cout << w << "\n";

	w = v;

	v[0] = 10;

	std::cout << w << "\n";
	std::cout << v << "\n";
	std::cout << v.dim() << "\n";

	Vector<float> z(v);
	z = z + w;
	z += w;
	std::cout << z << "\n";
	z = z - w;
	z -= w;
	std::cout << z << "\n";

	z = z * 5.f;
	z = 5.f * z;
	std::cout << z << "\n";

	z = z / 5.f;
	z = z / 5.f;
	std::cout << z << "\n";

	std::cout << -z << "\n";
	std::cout << +z << "\n";

	std::cout << (v * w) << "\n";
	std::cout << v.length() << "\n";

	Vector<float> a = { 6, -2 };
	Vector<float> b = { -4, 4 };
	std::cout << (a + b) << "\n";
	std::cout << (b + a) << "\n";

	Vector<float> c = { 2, 5 };
	Vector<float> d = { 7, 1 };
	std::cout << (c * d) << "\n";
	std::cout << c.length() << "\n";

	Vector<float> e = { 2, 2 };
	Vector<float> f = { 2, 0 };
	std::cout << e.angle(f) << "\n";

	Matrix<Vector<float>> M;

	M.addColumnVector(a);
	M.addColumnVector(b);

	std::cout << M << "\n";

	Matrix<Vector<float>> N;

	N.addColumnVector( Vector<float>{3, 6, 6} );
	N.addColumnVector( Vector<float>{1, 3, 8} );
	N.addColumnVector( Vector<float>{8, 0, 8} );

	std::cout << N << "\n";
	std::cout << N.getRowVector(2) << "\n";
	std::cout << N.getTranspose() << "\n";

	std::cout << N.getSubMatrix(0, 1) << "\n";

	std::cout << N.det() << "\n";

	return 0;
}
