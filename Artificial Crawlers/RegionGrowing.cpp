#include "RegionGrowing.h"

#include <queue>
int count = 0;

RegionGrowing::RegionGrowing( ITKImage3D::Pointer binaryImage )
{
    _image = binaryImage;
}



void RegionGrowing::execute( int radius )
{
    //createImageTest();
    if ( _image )
    {
		ITKIterator3D it(_image, _image->GetLargestPossibleRegion());

		for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
		{
			if (it.Get() == 255)
			{
				std::vector< Point3D >* group = new std::vector< Point3D >();

				ITKIterator3D::IndexType index = it.GetIndex();
				
				Point3D point(index[0], index[1], index[2]);

				grow(point, group, radius);

				_regions.push_back( group );
			}
        }
    }
}



std::vector< std::vector< Point3D >* > RegionGrowing::getRegions()
{
    return _regions;
}



void RegionGrowing::grow(Point3D point, std::vector< Point3D >* group, int radius)
{
	std::queue< Point3D > pointQueue;
	pointQueue.push(point);
	ITKImage3D::IndexType index;

	while (!pointQueue.empty())
	{
		point = pointQueue.front();
		pointQueue.pop();
		index[0] = point.x;
		index[1] = point.y;
		index[2] = point.z;

		if (_image->GetPixel(index) != 255)
			continue;

		_image->SetPixel(index, 0);

		group->push_back(point);

		int width, height, depth;
		ITKImage3D::SizeType size = _image->GetLargestPossibleRegion().GetSize();
		width = size[0];
		height = size[1];
		depth = size[2];

		int xStart = point.x >= radius ? point.x - radius : 0;
		int xEnd = point.x < ( width - radius ) ? point.x + radius : width - 1;

		int yStart = point.y >= radius ? point.y - radius : 0;
		int yEnd = point.y < (height - radius ) ? point.y + radius : height - 1;

		int zStart = point.z >= radius ? point.z - 1 : 0;
		int zEnd = point.z < (depth - radius ) ? point.z + radius : depth - 1;

		for (int z = zStart; z <= zEnd; z++)
		{
			for (int y = yStart; y <= yEnd; y++)
			{
				for (int x = xStart; x <= xEnd; x++)
				{
					ITKImage3D::IndexType otherIndex = { x, y, z };

					if (_image->GetPixel(otherIndex) == 255)
					{
						Point3D neighbor(x, y, z);
						pointQueue.push(neighbor);
					}
				}
			}
		}
	}
}

//void RegionGrowing::grow( Point3D point, Fault* fault )
//{
//    printf( "count: %d", ++count );
//    int width, height, depth;
//    _image->getSize( width, height, depth );
//
//    int xStart = point.x > 0 ? point.x - 1 : 0;
//    int xEnd   = point.x < width - 1 ? point.x + 1 : width - 1;
//    
//    int yStart = point.y > 0 ? point.y - 1 : 0;
//    int yEnd   = point.y < height - 1 ? point.y + 1 : height - 1;
//    
//    int zStart = point.z > 0 ? point.z - 1 : 0;
//    int zEnd   = point.z < depth - 1 ? point.z + 1 : depth - 1;
//
//    for ( int z = zStart; z <= zEnd; z++ )
//    {
//        for ( int y = yStart; y <= yEnd; y++ )
//        {
//            for ( int x = xStart; x <= xEnd; x++ )
//            {
//                if ( _image->getVoxel( x, y, z ) == 1 )
//                {
//                    _image->setVoxel( x, y, z, 0 );
//                    
//                    Point3D neighbor( x, y, z );
//                    fault->addPoint( neighbor );
//
//                    grow( neighbor, fault );
//                }
//            }
//        }
//    }
//    printf( "count: %d", --count );
//}


//void RegionGrowing::createImageTest()
//{
//    _image = new PFMImage( 60, 40 );
//
//    for ( int y = 7; y <= 15; y++ )
//        for ( int x = 35; x <= 45; x++ )
//            _image->setPixel( x, y, 1 );
//
//    for ( int y = 11; y <= 14; y++ )
//        for ( int x = 11; x <= 14; x++ )
//            _image->setPixel( x, y, 1 );
//
//    for ( int y = 10; y <= 15; y++ )
//        for ( int x = 17; x <= 18; x++ )
//            _image->setPixel( x, y, 1 );
//
//    _image->setPixel( 15, 13, 1 );
//    _image->setPixel( 16, 13, 1 );
//
//    for ( int y = 13; y <= 21; y++ )
//        for ( int x = 0; x <= 1; x++ )
//            _image->setPixel( x, y, 1 );
//
//    for ( int y = 31; y <= 35; y++ )
//        for ( int x = 21; x <= 48; x++ )
//            _image->setPixel( x, y, 1 );
//}