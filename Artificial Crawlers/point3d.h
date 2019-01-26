#ifndef POINT3D_H
#define POINT3D_H


class Point3D
{
public:
	double x;
	double y;
	double z;

	Point3D() : x(0), y(0), z(0) {}

	Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
/*
	Point3D(Point3D& p)
	{
		x = p.x;
		y = p.y;
		z = p.z;		
	}

*/
	bool operator ==(const Point3D &p)
	{
		if( x == p.x && y == p.y && z == p.z)
			return true;
		false;
	}

	Point3D& operator =(const Point3D &v)
    {
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
    }

};

#endif