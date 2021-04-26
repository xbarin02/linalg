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

	Matrix<float> M;

	M.addColumnVector(a);
	M.addColumnVector(b);

	std::cout << M << "\n";

	Matrix<float> N;

	N.addColumnVector( Vector<float>{3, 6, 6} );
	N.addColumnVector( Vector<float>{1, 3, 8} );
	N.addColumnVector( Vector<float>{8, 0, 8} );

	std::cout << N << "\n";
	std::cout << N.getRowVector(2) << "\n";
	std::cout << N.getTranspose() << "\n";

	std::cout << N.getSubMatrix(0, 1) << "\n";

	std::cout << N.det() << "\n";

	Vector<float> v1 = { 6, -2, 4 };
	Vector<float> v2 = { 4, 1, 5 };
	std::cout << (v1 * v2) << "\n";

	Matrix<float> A;
	A.addColumnVector(Vector<float>{1, 4, 7});
	A.addColumnVector(Vector<float>{2, 5, 8});
	A.addColumnVector(Vector<float>{3, 6, 9});

	Matrix<float> B;
	B.addColumnVector(Vector<float>{4, 1, 3});
	B.addColumnVector(Vector<float>{2, 7, 2});

	std::cout << (A * B) << "\n";

	Matrix<float> C;
	C.addColumnVector(Vector<float>{7, 5, 3});
	C.addColumnVector(Vector<float>{6, 4, 2});

	std::cout << (B + C) << "\n";

	std::cout << Matrix<float>::identity(4) << "\n";

	std::cout << A << "\n";
	A.swapRowVectors(0, 2);
	A.addVectorToRow(1, Vector<float>{2, 1, 0});
	A.subVectorFromRow(1, Vector<float>{1, 1, 1});
	std::cout << A << "\n";

	Matrix<float> D;
	D.addColumnVector(Vector<float>{1, 4, 7});
	D.addColumnVector(Vector<float>{2, 5, 8});
	D.addColumnVector(Vector<float>{3, 6, 9});
	std::cout << D << "\n";
	std::cout << D.det() << "\n";
	D.toRref();
	std::cout << D << "\n";

	Matrix<float> E;
	E.addColumnVector(Vector<float>{0, 0, 0, 0});
	E.addColumnVector(Vector<float>{0, 0, 0, 0});
	E.addColumnVector(Vector<float>{5, 6, 7, 5});
	E.addColumnVector(Vector<float>{5, 6, 8, 5});
	E = E.getTranspose();
	std::cout << E << "\n";
	std::cout << E.getBasis() << "\n";
	std::cout << E.det() << "\n";
	E.toRref();
	std::cout << E << "\n";

	std::cout << E.getColumnVector(0).is_pivot() << " " << E.getColumnVector(0).getPivotEntry() << "\n";
	std::cout << E.getColumnVector(1).is_pivot() << " " << E.getColumnVector(1).getPivotEntry() << "\n";
	std::cout << E.getColumnVector(2).is_pivot() << " " << E.getColumnVector(2).getPivotEntry() << "\n";
	std::cout << E.getColumnVector(3).is_pivot() << " " << E.getColumnVector(3).getPivotEntry() << "\n";
	std::cout << "here" << "\n\n";

	// rref(F) = I
	Matrix<float> F;
	F.addColumnVector(Vector<float>{1, 4, 7});
	F.addColumnVector(Vector<float>{4, 2, 8});
	F.addColumnVector(Vector<float>{3, 1, 1});
	std::cout << F << "\n";
	std::cout << F.det() << "\n";
// 	F.toRref();
// 	std::cout << F << "\n";
	Matrix<float> G = F.getInverse();
	std::cout << G << "\n";
	Matrix<float> H = F * G;
	H.round();
	std::cout << H << "\n";
	std::cout << F.getRref() << "\n";


	return 0;
}
