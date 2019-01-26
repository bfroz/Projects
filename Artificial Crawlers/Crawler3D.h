#include <itkImage.h>

#ifndef Crawler3D_H
#define Crawler3D_H
typedef short PixelType;
typedef itk::Image<PixelType, 3> ITKImage3D;

class Crawler3D{
public:
	int energy;
	int min_energy;
	int time;
	std::vector<float> move_hist;
	std::vector<ITKImage3D::IndexType> move_hist_index;
	ITKImage3D::IndexType position;
	Crawler3D(int l, int ml, ITKImage3D::IndexType p){
		energy = l, min_energy = ml, position = p;
	}
	Crawler3D(){}
	void move(ITKImage3D::IndexType p){
		attenuate();
		float x = position.GetIndex()[0];
		float y = position.GetIndex()[1];
		float newX = p.GetIndex()[0];
		float newY = p.GetIndex()[0];
		move_hist.push_back(atan2(y - newY, x - newX));
		move_hist_index.push_back(p);
		position = p;
	}
	void setLife(int i){
		energy = i;
	}
	void attenuate(){
		energy -= 1;
	}
	void accentuate(int hp){
		energy += hp;
	}
	std::vector<float> getMoveHistory(){
		return move_hist;
	}
	std::vector<ITKImage3D::IndexType> getMoveHistoryIndex(){
		return move_hist_index;
	}
	std::vector<ITKImage3D::IndexType> getMoveHistoryIndexWithPosition(){
		if (move_hist_index.size() == 0){
			move_hist_index.push_back(position);
		}
		return move_hist_index;
	}
};
#endif