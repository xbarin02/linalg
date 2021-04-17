#include "Vector.h"

#include <ostream>

using namespace LinAlg;

int main()
{
	Vector<float> v(5);

	std::cout << v << "\n";

	v[2] = 3.14;
	std::cout << v[2] << "\n";

	std::cout << v << "\n";

	Vector<float> w(v);

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

	return 0;
}
