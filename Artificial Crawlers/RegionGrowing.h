#ifndef REGIONGROWING_H
#define REGIONGROWING_H

#include <vector>

#include "DEFINICOES.h"
#include "point3d.h"

class RegionGrowing
{
public:
    RegionGrowing( ITKImage3D::Pointer binaryImage );

    void execute( int radius = 1 );

    std::vector< std::vector< Point3D >* > getRegions();

private:
    //void grow( Point3D point, std::vector< Point3D >* group );
	void grow(Point3D point, std::vector< Point3D >* group, int radius = 1);

    void createImageTest();

    ITKImage3D::Pointer _image;
    std::vector< std::vector< Point3D >* > _regions;
};


#endif
