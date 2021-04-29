#ifndef FLOAT_H
#define FLOAT_H

#include <cmath>

namespace LinAlg {

template <typename T>
class Float
{
public:
	static bool nearly(T a, T b = 0, T epsilon = 1e-6)
	{
		float diff = std::fabs(a - b);

		return diff < epsilon;
	}
};

}

#endif
