#include "Point.h"
#include <vector>
#include <numeric>
#include <functional>
#include <iostream>

typedef std::vector< PointN > Curve_T;

//The "double error" variable is a radius that gives a margin of how much the elements must be equal.
//In other words: if the error is 0 (zero), the elements must be equal; If it's 1.5, the elements must be "almost" equals, with a +-1.5 margin error;
Real_T point_distance(const PointN& a, const PointN& b)
{
	return hypot(a.x - b.x, a.y - b.y);
}

Real_T curve_distance(const Curve_T& c1, const Curve_T& c2)
{
	if (c1.size() != c2.size()) throw std::invalid_argument("size mismatch");
	return std::inner_product(c1.begin(), c1.end(), c2.begin(), Real_T(0), std::plus< Real_T >(), point_distance);
}
Real_T curve_jaccard_distance(const Curve_T& c1, const Curve_T& c2, double error){
	// d(c1,c2) = (|c1Uc2| - |c1Ac2|) / |c1Uc2|
	double U = c1.size(), A = 0.0;
	for (int i = 0; i < U; i++){
		if ((c1.at(i).y + error) >= c2.at(i).y && (c1.at(i).y - error) <= c2.at(i).y) A++;
	}
	return (U - A) / U;
}

Real_T curve_simple_matching(const Curve_T& c1, const Curve_T& c2, double error){
	// d(c1,c2) = a / |(a+b+c)
	Real_T a = 0, b = 0, c = 0;
	for (int i = 0; i < c1.size(); i++){
		if ((c1.at(i).y + error) >= c2.at(i).y && (c1.at(i).y - error) <= c2.at(i).y) a++;
		else{
			b++; c++;
		}
	}
	return (Real_T)1 - (a / (a + b + c));
}

Real_T curve_hamming_distance(const Curve_T& c1, const Curve_T& c2, double error){
	int equal_bits = 0;
	for (int i = 0; i < c1.size(); i++){
		if ((c1.at(i).y + error) >= c2.at(i).y && (c1.at(i).y - error) <= c2.at(i).y) equal_bits++;
	}
	return equal_bits;
}
//Chebyshev distance : D_{\rm Chebyshev}(p,q) := \max_i(|p_i - q_i|).\ 
Real_T curve_chebyshev_distance(const Curve_T& c1, const Curve_T& c2, double error){
	Real_T cumulator = 0;
	Real_T max = 0;
	for (int i = 0; i < c1.size(); i++){
		cumulator = sqrt(pow((c1.at(i).x - c2.at(i).x), 2) + pow((c1.at(i).y - c2.at(i).y), 2));
		if (cumulator > max){
			max = cumulator;
		}
	}
	return max;
}
//Manhattan distance : d_1(\mathbf{p}, \mathbf{q}) = \|\mathbf{p} - \mathbf{q}\|_1 = \sum_{i=1}^n |p_i-q_i|,
Real_T curve_manhattan_distance(const Curve_T& c1, const Curve_T& c2, double error){
	Real_T cumulator = 0;
	for (int i = 0; i < c1.size(); i++){
		cumulator += abs(c1.at(i).x - c2.at(i).x) + abs(c1.at(i).y - c2.at(i).y);
	}
	return cumulator;
}

class Curve{
public: Curve_T curve;

public: Curve(){}

		void setElement(float x, float y){
			curve.push_back(PointN(x, y));
		}
		void printElements(){
			normalize(&curve);
			for (int i = 0; i < curve.size(); i++){
				std::cout << "X: " << curve.at(i).x << " Y: " << curve.at(i).y << "\n";
			}
		}
		float getArea(){
			return polygonArea();
		}
		Curve_T getCurve(){
			normalize(&curve);
			return curve;
		}
private:
	float polygonArea()
	{
		normalize(&curve);
		curve.push_back(PointN(curve.size(), 0));
		curve.push_back(PointN(0, 0));

		float area = 0;         // Accumulates area in the loop
		int j = curve.size() - 1;  // The last vertex is the 'previous' one to the first

		for (int i = 0; i < curve.size(); i++)
		{
			area = area + (curve.at(j).x + curve.at(i).x) * (curve.at(j).y - curve.at(i).y);
			j = i;  //j is previous vertex to i
		}
		return area / 2;
	}

	void normalize(Curve_T * vetor){
		float new_y, lower, higher;
		lower = getLowerY(*vetor);
		higher = getHigherY(*vetor);
		for (int i = 0; i < vetor->size(); i++){
			if (higher == lower) new_y = 0;
			else new_y = (vetor->at(i).y - lower) / (higher - lower);
			vetor->at(i).y = new_y;
		}
	}
	float getLowerY(Curve_T  vetor){
		float lower = vetor.at(0).y;
		for (int i = 0; i < vetor.size(); i++){
			if (vetor.at(i).y < lower) lower = vetor.at(i).y;
		}
		return lower;
	}
	float getHigherY(Curve_T  vetor){
		float higher = vetor.at(0).y;
		for (int i = 0; i < vetor.size(); i++){
			if (vetor.at(i).y > higher) higher = vetor.at(i).y;
		}
		return higher;
	}

};