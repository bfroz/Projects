#include <initializer_list>

typedef double Real_T;

class PointN
{
public:
	PointN(){}
	PointN(std::initializer_list<Real_T> args) :x(args.begin()[0]), y(args.begin()[1]){}
	PointN(const Real_T& xx, const Real_T& yy) :x(xx), y(yy){}
	Real_T x, y;
};