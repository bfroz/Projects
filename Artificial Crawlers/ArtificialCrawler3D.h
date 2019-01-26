#include "Crawler3D.h"
//#include "Curve.h"
#include "RegionGrowing.h"

//#include <itkSobelOperator.h>
#include <ctime>
#include "itkRegionOfInterestImageFilter.h"
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include "Util.h"

#ifndef ArtificialCrawler3D_H
#define ArtificialCrawler3D_H
typedef itk::RegionOfInterestImageFilter< ITKImage3D, ITKImage3D > RegionFilterType2;
typedef itk::NeighborhoodIterator<ITKImage3D> NeighborhoodIterator2;
typedef itk::NeighborhoodInnerProduct<ITKImage3D> NeighborhoodInnerProduct2;
//typedef itk::NeighborhoodConnectedImageFilter<ITKImage3D,ITKImage3D > ConnectedFilterType;

//RULES
//unsigned const rule1 = 1;
//unsigned const rule2 = 2;
//unsigned const rule3 = 3;

class ArtificialCrawler3D{
	ITKImage3D::Pointer environment;
	int min_environment_value;
	int sizex, sizey, sizez;
	std::vector<Crawler3D>* population;
	int start_energy;
	Curve curveNumberCrawlers;
	Curve curveColony;
	Curve curveScaleDistribuition;
	Curve curveHabitantsSettlement;
	std::vector<std::vector<Crawler3D>> global_clusters;
	int cont_cycles;
public:
	ArtificialCrawler3D(){}
	//
	ArtificialCrawler3D(ITKImage3D::Pointer env, int s_energy){
		environment = env;
		sizex = env->GetLargestPossibleRegion().GetSize()[0]; 
		sizey = env->GetLargestPossibleRegion().GetSize()[1];
		sizez = env->GetLargestPossibleRegion().GetSize()[2];
		start_energy = s_energy;
		population = new std::vector<Crawler3D>();
		min_environment_value = getMinimumValue3D(); //finds the minimum value - image background value.
		createPopulation(); //Creates the population.
		cont_cycles = 0;
	}
	bool sortIndexByPixelIntensity(ITKImage3D::IndexType i, ITKImage3D::IndexType j) { return environment->GetPixel(i) > environment->GetPixel(j); }

	//createPopulation() - creates the initial population based on pixel value. Each element in the vector "population" is a Crawler object.
	//The condition to the creation of a Crawler is the minimum value of the image - a Crawler will be created on a non-minimum pixel value position of the environment.
	void createPopulation(){
		itk::ImageRegionIterator<ITKImage3D> itE(environment, environment->GetRequestedRegion());
		while (!itE.IsAtEnd())
		{
			if (itE.Get() > min_environment_value){
				Crawler3D crawler = Crawler3D(start_energy, 0, itE.GetIndex());
				population->push_back(crawler);
			}
			++itE;
		}
	}

	//getBinaryPopulation() - returns an 3D binary image showing all Crawlers position - all white pixels are Crawlers from the actual population.
	ITKImage3D::Pointer getBinaryPopulation(){
		ITKImage3D::Pointer output = ITKImage3D::New();
		ITKImage3D::SizeType size;
		size[0] = sizex;
		size[1] = sizey;
		size[2] = sizez;
		ITKImage3D::IndexType start;
		start.Fill(0);
		ITKImage3D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		output->FillBuffer(0);
		ITKImage3D::IndexType index;
		index[0] = -1;
		index[1] = -1;
		index[2] = -1;
		for (unsigned i = 0; i < population->size(); i++){
			//output->SetPixel(index, 255);	
			index = population->at(i).position;
			if (index[0] < sizex && index[1] < sizey)
				output->SetPixel(index, 255);
			//std::cout << "\n" << i;
		}
		return output;
	}
	
	//getPopulation() - returns an 3D image showing all Crawlers position - all non-minimum value pixels are Crawlers from the actual population.
	ITKImage3D::Pointer getPopulation(){
		ITKImage3D::Pointer output = ITKImage3D::New();
		ITKImage3D::SizeType size;
		size[0] = sizex;
		size[1] = sizey;
		size[2] = sizez;
		ITKImage3D::IndexType start;
		start.Fill(0);
		ITKImage3D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		output->FillBuffer(min_environment_value);
		ITKImage3D::IndexType index;
		index[0] = -1;
		index[1] = -1;
		index[2] = -1;
		for (unsigned i = 0; i < population->size(); i++){
			index = population->at(i).position;
			if (index[0] < sizex && index[1] < sizey)
				output->SetPixel(index, environment->GetPixel(index));
		}
		return output;
	}
	
	//start(int cycles) - Starts running Crawlers. It runs until 'cycles' is reached; I.e. If 'cycles' = 5, it runs the first 5 cycles of the ACrawlers
	//The rules are descrived bellow.
	//Artificial Crawlers Life
	//Rules
	//- If the maximum intensity of one of the 8-neighborhood is lower than the current Crawler, it keeps unchanged (1)
	//- If the maximum intensity of one of the 8-neighborhood is higher than the current Crawler, it moves to the higher place (2)
	//- If there is at least one place on the 8-neighborhood higher than the current Crawler, moves to the one that another Crawler moved before (3)
	//-- If another Crawler is on the same place, let them 'fight'
	//The energy of each Crawler reduces when it moves
	//The amount of absorbed energy from another Crawler is a fixed amount (20%)
	//The amount of absorbed energy from environment is a fixed amount (10%)
	void start(int cycles){
		std::vector<ITKImage3D::IndexType> neighborhood;
		std::vector<ITKImage3D::IndexType> higher_pixels;
		int dead_cont = 0;
		int rule = 0;
		cont_cycles = cycles;
		int maximum_neighborhood_intensiy; //The maximum of intensity of the 8 neighborhood
		bool hasCrawler = false;
		Crawler3D enemy;
		int enemy_position;
		int randPosition;
		//Begin the Cycles
		for (unsigned int x = 0; x < cycles; x++){
			curveNumberCrawlers.setElement(x, population->size());
			dead_cont = 0;
			if (population){
				for (int i = population->size() - 1; i > 1; i--){
					if (population->at(i).energy <= 0) continue; //If the actual Crawler has no energy, do not move
					//Get the voxel neighborhood of an Index
					neighborhood = getVoxelNeighborhood(population->at(i).position);
					for (int j = 0; j < neighborhood.size(); j++){
						if (environment->GetPixel(neighborhood.at(j)) > environment->GetPixel(population->at(i).position)){
							//Get all higher pixels on neighborhood
							higher_pixels.push_back(neighborhood.at(j));
						}
					}
					//Compare the higher pixels to find out which rule to use
					if (higher_pixels.size() == 0) rule = rule1;
					else if (higher_pixels.size() == 1) rule = rule2;
					else rule = rule3;

					switch (rule)
					{
					case 1:
						//Just...don't move.
						break;
					case 2:
						//Move to the only one higher.
						hasCrawler = hasCrawlerSet(higher_pixels.at(0), &enemy, &enemy_position);
						if (hasCrawler){
							if (population->at(i).energy >= enemy.energy){
								//If energy of the actual Crawler is higher or equal than his actual enemy
								population->at(i).accentuate((int)enemy.energy * 0.2);
								population->at(i).move(higher_pixels.at(0));
								population->erase(population->begin() + enemy_position);
							}
							else{
								//If energy of the enemy is higher
								population->at(enemy_position).accentuate((int)population->at(i).energy * 0.2);
								population->erase(population->begin() + i);
							}
							dead_cont++;
						}
						else population->at(i).move(higher_pixels.at(0));
						break;
					case 3:
						//Move to the one which has a Crawler marked. If there is no Crawler, move to the first one.
						//First, sort the neighborhood from the higher to the lower pixel
						//std::sort(higher_pixels.begin(), higher_pixels.end(), sortIndexByPixelIntensity); //TODO: Corrigir isso
						for (int k = 0; k < higher_pixels.size(); k++){
							//Check if exists another Crawler in the higher pixel position
							hasCrawler = hasCrawlerSet(higher_pixels.at(k), &enemy, &enemy_position);
							if (hasCrawler){
								if (k == higher_pixels.size() - 1){
									//If it's the last higher Crawler, sort a position from k=0 to k=higher_pixels.size()-1
									if (population->at(i).energy >= enemy.energy){
										//If energy of the actual Crawler is higher or equal than his actual enemy
										population->at(i).accentuate((int)enemy.energy * 0.2);
										population->at(i).move(higher_pixels.at(k));
										population->at(i).accentuate(std::abs(environment->GetPixel(higher_pixels.at(0))) * 0.1); //Absorbs 10% of the environment energy //TODO: parameterize this to absorbs ?%.
										population->erase(population->begin() + enemy_position);
									}
									else{
										//If energy of the enemy is higher
										population->at(enemy_position).accentuate((int)population->at(i).energy * 0.2);
										population->erase(population->begin() + i);
									}
									dead_cont++;
									break;
								}
								else continue;
							}
							else
							{
								//Else, there is no Crawler in the current higher pixel
								population->at(i).move(higher_pixels.at(k));
								population->at(i).accentuate(std::abs(environment->GetPixel(higher_pixels.at(0))) * 0.1); //Absorbs 10% of the environment energy //TODO: parameterize this to absorbs ?%.

								break;
							}

						}
					default:
						break;
					}
					neighborhood.clear();
					higher_pixels.clear();
				}
				curveHabitantsSettlement.setElement(x, dead_cont);
			}
			eliminateDeadCrawlers(); //TODO: Sets a minimum threshold //Fixed in 0
		}
		global_clusters = ConnectedElementsByRadius(1);
		//Get Characteristics
		createCurveScaleDistribuition(10);
		createCurveColony(10);
	}

	//start() - Starts running Crawlers. It runs until all ACrawlers use rule 1 (none of them moved)
	//The rules are descrived bellow.
	//Artificial Crawlers Life
	//Rules
	//- If the maximum intensity of one of the 8-neighborhood is lower than the current Crawler, it keeps unchanged (1)
	//- If the maximum intensity of one of the 8-neighborhood is higher than the current Crawler, it moves to the higher place (2)
	//- If there is at least one place on the 8-neighborhood higher than the current Crawler, moves to the one that another Crawler moved before (3)
	//-- If another Crawler is on the same place, let them 'fight'
	//The energy of each Crawler reduces when it moves
	//The amount of absorbed energy from another Crawler is a fixed amount (20%)
	//The amount of absorbed energy from environment is a fixed amount (10%)
	void start(){
		std::vector<ITKImage3D::IndexType> neighborhood;
		std::vector<ITKImage3D::IndexType> higher_pixels;
		int dead_cont = 0;
		int rule = 0;
		int maximum_neighborhood_intensiy; //The maximum of intensity of the 8 neighborhood
		bool hasCrawler = false;
		Crawler3D enemy;
		int enemy_position;
		int randPosition;
		int cont_moves;
		
		//Begin the Cycles
		while (true){			
			cont_moves = 0;
			curveNumberCrawlers.setElement(cont_cycles, population->size());
			dead_cont = 0;
			if (population){
				for (int i = population->size() - 1; i > 1; i--){
					if (population->at(i).energy <= 0) continue; //If the actual Crawler has no energy, do not move
					//Get the voxel neighborhood of an Index
					neighborhood = getVoxelNeighborhood(population->at(i).position);
					for (int j = 0; j < neighborhood.size(); j++){
						if (environment->GetPixel(neighborhood.at(j)) > environment->GetPixel(population->at(i).position)){
							//Get all higher voxels on neighborhood
							higher_pixels.push_back(neighborhood.at(j));
						}
					}
					//Compare the higher voxels to find out which rule to use
					if (higher_pixels.size() == 0) rule = rule1;
					else if (higher_pixels.size() == 1) rule = rule2;
					else rule = rule3;

					switch (rule)
					{
					case 1:
						//Just...don't move.
						break;
					case 2:
						cont_moves++;
						//Move to the only one higher.
						hasCrawler = hasCrawlerSet(higher_pixels.at(0), &enemy, &enemy_position);
						if (hasCrawler){
							if (population->at(i).energy >= enemy.energy){
								//If energy of the actual Crawler is higher or equal than his actual enemy
								population->at(i).accentuate((int)enemy.energy * 0.2); //Absorbs a little energy from enemy //TODO: absorbs all enemy's energy
								population->at(i).move(higher_pixels.at(0));
								population->erase(population->begin() + enemy_position);
								population->at(i).accentuate(std::abs(environment->GetPixel(higher_pixels.at(0))) * 0.1); //Absorbs 10% of the environment energy //TODO: parameterize this to absorbs ?%.

							}
							else{
								//If energy of the enemy is higher
								population->at(enemy_position).accentuate((int)population->at(i).energy * 0.2);
								population->erase(population->begin() + i);

							}
							dead_cont++;
						}
						else {
							population->at(i).move(higher_pixels.at(0));
							population->at(i).accentuate(std::abs(environment->GetPixel(higher_pixels.at(0))) * 0.1); //Absorbs 10% of the environment energy //TODO: parameterize this to absorbs ?%.
						}
						break;
					case 3:
						cont_moves++;
						//Move to the one which has a Crawler marked. If there is no Crawler, move to the first one.
						//First, sort the neighborhood from the higher to the lower pixel
						//std::sort(higher_pixels.begin(), higher_pixels.end(), sortIndexByPixelIntensity); //TODO: Do this correctly (optional)
						for (int k = 0; k < higher_pixels.size(); k++){
							//Check if exists another Crawler in the higher pixel position
							hasCrawler = hasCrawlerSet(higher_pixels.at(k), &enemy, &enemy_position);
							if (hasCrawler){
								if (k == higher_pixels.size() - 1){
									//If it's the last higher Crawler, sort a position from k=0 to k=higher_pixels.size()-1
									if (population->at(i).energy >= enemy.energy){
										//If energy of the actual Crawler is higher or equal than his actual enemy
										population->at(i).accentuate((int)enemy.energy * 0.2);
										population->at(i).move(higher_pixels.at(k));
										population->erase(population->begin() + enemy_position);
									}
									else{
										//If energy of the enemy is higher
										population->at(enemy_position).accentuate((int)population->at(i).energy * 0.2);
										population->erase(population->begin() + i);
									}
									dead_cont++;
									break;
								}
								else continue;
							}
							else
							{
								//Else, there is no Crawler in the current higher pixel
								population->at(i).move(higher_pixels.at(k));
								break;
							}

						}
					default:
						break;
					}
					neighborhood.clear();
					higher_pixels.clear();
				}
				/*saveImage(getPopulation(), ss.str().c_str());*/
				curveHabitantsSettlement.setElement(cont_cycles, dead_cont); //Sets the habitantsSettlement elements - Number of dead crawlers per cycle.
			}
			if (cont_moves == 0) break;
			cont_cycles++;
			eliminateDeadCrawlers(); //TODO: Sets a minimum threshold //Fixed in 0
		}
		//Create global
		global_clusters = ConnectedElementsByRadius(1); //It defines the clusters/colony with 1 radius range; 
		// - These colonies are all elements connected in the end of ACrawlers cycles
		//Get Characteristics
		createCurveScaleDistribuition(10); //Creates the first 10 cycles curve //TODO: use some parameterization, global or local number to define the X size of the curve
		createCurveColony(10); //Creates the first 10 cycles curve //TODO: use some parameterization, global or local number to define the X size of the curve
	}
	//getCyclesCont() - returns the number of cycles (only after 'start(...)' function ends)
	int getCyclesCont(){
		return cont_cycles;
	}
	//getPopulationSize() - returns the actual population size.
	int getPopulationSize(){
		return population->size();
	}

	
	//QTClustering(...) - Apply the QT Clustering in the population.
	//--diameter_size = maximum diameter to be used in the Cluster formation - Higher the value, more elements inside the Clusters.
	//--minimum_cardinality = the minimum number of elements inside each Cluster.
	//--minimum_quality = TODO: implement minimum quality.
	std::vector<std::vector<Crawler3D>> QTClustering(float diameter_size, int minimum_cardinality, int minimum_quality){
		std::vector<std::vector<Crawler3D>> ClustersResult = std::vector<std::vector<Crawler3D>>();
		std::vector<Crawler3D> aux_population = *population;

		QTClusteringRecursiveCall(&aux_population, &ClustersResult, diameter_size, minimum_quality);
		//Check cardinality. Erase Clusters with cardinality lower than minimum_cardinality
		for (int i = ClustersResult.size() - 1; i >= 0; i--){
			if (ClustersResult.at(i).size() < minimum_cardinality){
				ClustersResult.erase(ClustersResult.begin() + i);
			}
		}
		return ClustersResult;
	}
	
	//ConnectedElementsByRadius(...) - returns all colonies within a maximum 'radius'.
	std::vector<std::vector<Crawler3D>> ConnectedElementsByRadius(int radius){
		std::vector<std::vector<Crawler3D>> output;
		std::vector<std::vector<Point3D>*> outroutput;
		ITKImage3D::Pointer binary_population = getBinaryPopulation(); //Gets the binary image representation of the population
		RegionGrowing rg = RegionGrowing(binary_population);
		rg.execute(radius);
		outroutput = rg.getRegions();
		output = ConvertGroupOfPointsToColony(outroutput);
		return output;
	}

	void printCurveNumberOfCrawlers(){
		curveNumberCrawlers.printElements();
	}
	void printCurveColony(int size){
		createCurveColony(size);
		curveColony.printElements();
	}
	void printCurveScaleDistribution(int size){
		createCurveScaleDistribuition(size);
		curveScaleDistribuition.printElements();
		std::cout << "Area: " << curveScaleDistribuition.getArea();
	}
	void printCurveHabitantsSettlement(){
		curveHabitantsSettlement.printElements();
	}
	float getCurveScaleDistributionArea(){
		return curveScaleDistribuition.getArea();
	}
	float getCurveColonyFormationnArea(){
		return curveColony.getArea();
	}
	float getCurveNumberCrawlersArea(){
		return curveNumberCrawlers.getArea();
	}
	float getCurveHabitantsSettlementArea(){
		return curveHabitantsSettlement.getArea();
	}
	Curve_T getCurveScaleDistribution(){
		/*std::cout << "\n Scale: \n";
		curveScaleDistribuition.printElements();*/
		return curveScaleDistribuition.getCurve();
	}
	Curve_T getCurveColonyFormation(){
		/*std::cout << "\n Colony: \n";
		curveColony.printElements();*/
		return curveColony.getCurve();
	}
	Curve_T getCurveNumberCrawlers(){
		/*std::cout << "\n Number: \n";
		curveNumberCrawlers.printElements();*/
		return curveNumberCrawlers.getCurve();
	}
	Curve_T getCurveHabitantsSettlement(){
		/*std::cout << "\n Habitants: \n";
		curveHabitantsSettlement.printElements();*/
		return curveHabitantsSettlement.getCurve();
	}

private:
	//getMinimumValue3D() - get the background value of the environment - the minimum voxel value.
	int getMinimumValue3D(){
		typedef itk::MinimumMaximumImageCalculator <ITKImage3D> ImageCalculatorFilterType;
		ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
		imageCalculatorFilter->SetImage(environment);
		imageCalculatorFilter->ComputeMinimum();
		return imageCalculatorFilter->GetMinimum();
	}

	//hasCrawlerSet(...) - returns true if has a Crawler in the "index" parameter. If has Crawler, "enemy" and "enemy_position" will be set.
	bool hasCrawlerSet(ITKImage3D::IndexType index, Crawler3D* enemy, int * enemy_position){
		for (unsigned int i = 0; i < population->size(); i++){
			if (population->at(i).position == index){
				*enemy = population->at(i);
				*enemy_position = i;
				return true;
			}
		}
		return false;
	}

	//eliminateDeadCrawlers() - eleminates Crawlers with minimum energy. //TODO: create a minimum threshold to represents minimum energy.
	void eliminateDeadCrawlers(){
		for (unsigned int i = population->size() - 1; i < 0; i--){
			if (population->at(i).energy <= 0){
				population->erase(population->begin() + i);
			}
		}
	}
	
	//getVoxelNeighborhood(...) - get the nearest elements of the "crawler_index".
	//It returns a vector with all neighborhood elements.
	std::vector<ITKImage3D::IndexType> getVoxelNeighborhood(ITKImage3D::IndexType crawler_index){
		std::vector<ITKImage3D::IndexType> neighborhood;
		RegionFilterType2::Pointer Regionfilter = RegionFilterType2::New();
		ITKImage3D::RegionType desiredRegion;
		ITKImage3D::Pointer neighborhoodImage = ITKImage3D::New();
		itk::Size<3> size;
		itk::Index<3> index;
		itk::Index<3> new_index;
		if (crawler_index[0] > 0){
			index[0] = crawler_index[0] - 1;
			size[0] = 3;
		}
		else{
			index[0] = crawler_index[0];
			size[0] = 2;
		}
		if (crawler_index[1] > 0){
			index[1] = crawler_index[1] - 1;
			size[1] = 3;
		}
		else{
			index[1] = crawler_index[1];
			size[1] = 2;
		}
		if (crawler_index[2] > 0){
			index[2] = crawler_index[2] - 1;
			size[2] = 3;
		}
		else{
			index[2] = crawler_index[2];
			size[2] = 2;
		}
		if ((index[0] + size[0]) >= sizex){
			size[0] = sizex - index[0];
		}
		if ((index[1] + size[1]) >= sizey){
			size[1] = sizey - index[1];
		}
		if ((index[2] + size[2]) >= sizez){
			size[2] = sizez - index[2];
		}
		desiredRegion.SetSize(size);
		desiredRegion.SetIndex(index);
		Regionfilter->SetRegionOfInterest(desiredRegion);
		Regionfilter->SetInput(environment);
		Regionfilter->Update();
		neighborhoodImage = Regionfilter->GetOutput();

		itk::ImageRegionIterator<ITKImage3D> it(neighborhoodImage, neighborhoodImage->GetRequestedRegion());
		while (!it.IsAtEnd())
		{
			new_index[0] = it.GetIndex()[0] + index[0];
			new_index[1] = it.GetIndex()[1] + index[1];
			new_index[2] = it.GetIndex()[2] + index[2];
			if (new_index != crawler_index){
				neighborhood.push_back(new_index);
			}
			++it;
		}
		return neighborhood;
	}
	
	//QTClusteringRecursiveCall(...) - Recursive calling of the QT Clustering.
	//"candidates" and "ClustersResult" pointers are updated inside the method.
	void QTClusteringRecursiveCall(std::vector<Crawler3D> *candidates, std::vector<std::vector<Crawler3D>> *ClustersResult, float diameter_size, int minimum_quality){
		std::vector<Crawler3D> CurrentCluster = std::vector<Crawler3D>(); //Current Cluster. This element will save the actual cluster of the QT iterations
		std::vector<std::vector<Crawler3D>> TempClustersCandidates = std::vector<std::vector<Crawler3D>>();
		float distance;
		if (candidates->size() > 0){ //If there is candidates, do QT
			for (unsigned i = 0; i < candidates->size(); i++){
				CurrentCluster.push_back(candidates->at(i)); //Put the first candidate in a cluster
				//Check every candidate. The candidates that have a distance from the current cluster lower or equal than diameter_size will be a new participant of the CurrentCluster
				for (unsigned j = 0; j < candidates->size(); j++){
					if (j == i) continue; //If its the CurrentCluster seed, continue, because it is already in the CurrentCluster
					distance = distancePointToGroup(CurrentCluster, candidates->at(j));
					//If the CurrentCluster distance to the actual candidate is equal or lower than the diameter_size, put the actual candidate in the cluster
					if (distance <= diameter_size){
						CurrentCluster.push_back(candidates->at(j));
					}
				}
				//Add the CurrentCluster in the vector of clusters
				TempClustersCandidates.push_back(CurrentCluster);
				CurrentCluster.clear();
			}
			//Get the maximum cardinality cluster in TempClustersCandidates TODO: Check the minimum_quality threshold to choose the best candidate
			if (TempClustersCandidates.size() > 0){
				std::vector<Crawler3D> higherCluster = CurrentCluster;
				for (unsigned k = 0; k < TempClustersCandidates.size(); k++){
					if (higherCluster.size() < TempClustersCandidates.at(k).size()){
						higherCluster = TempClustersCandidates.at(k);
					}
				}
				//Add the best Cluster in the ClustersResult
				ClustersResult->push_back(higherCluster);
				//Remove the candidates inside the best Cluster the from *candidates
				for (unsigned x = 0; x < higherCluster.size(); x++){
					int position = -1;
					for (unsigned z = 0; z < candidates->size(); z++){
						if (candidates->at(z).position == higherCluster.at(x).position){
							position = z;
							break;
						}
					}
					if (position >= 0){
						candidates->erase(candidates->begin() + position);
					}
				}
				//Call QTClusteringRecursiveCall(------)
				QTClusteringRecursiveCall(candidates, ClustersResult, diameter_size, minimum_quality);
			}

		}
	}

	//distancePointToGroup(...) - calculates the distance between a Crawler and a group of Crawler.
	float distancePointToGroup(std::vector<Crawler3D> group, Crawler3D point){
		float distance = 99999999999; //TODO: puts MAX INT
		float aux;
		for (unsigned int i = 0; i < group.size(); i++){
			aux = distancePointToPoint(group.at(i), point);
			if (aux < distance){
				distance = aux;
			}
		}
		return distance;
	}

	//distancePointToPoint(...) - calculates the distance between two Crawlers.
	float distancePointToPoint(Crawler3D p1, Crawler3D p2){
		float distance = 0;
		return sqrt(pow(((float)p2.position[0] - (float)p1.position[0]), 2) + pow(((float)p2.position[1] - (float)p1.position[1]), 2) + pow(((float)p2.position[2] - (float)p1.position[2]), 2));
	}

	//ConvertGroupOfPointsToColony(...) - Converts a vector of vector<Point3D> into a vector of vector<Crawler3D>
	std::vector<std::vector<Crawler3D>> ConvertGroupOfPointsToColony(std::vector<std::vector<Point3D>*> groups){
		std::vector<std::vector<Crawler3D>> output;
		std::vector<Crawler3D> colony;
		ITKImage3D::IndexType index;
		for (int i = 0; i < groups.size(); i++){
			for (int j = 0; j < groups.at(i)->size(); j++){
				index[0] = groups.at(i)->at(j).x;
				index[1] = groups.at(i)->at(j).y;
				index[2] = groups.at(i)->at(j).z;
				for (int k = 0; k < population->size(); k++){					
					if (index == population->at(k).position){
						colony.push_back(population->at(k));
						break;
					}
				}				
			}
			output.push_back(colony);
			colony.clear();
		}
		return output;
	}

	void saveImage(ITKImage3D::Pointer input, const char* path)
	{
		typedef  itk::ImageFileWriter<ITKImage3D> WriterType;
		WriterType::Pointer Writer = WriterType::New();
		std::stringstream ss;
		ss << path << ".mhd";
		Writer->SetFileName(ss.str().c_str());
		Writer->SetInput(input);
		Writer->Update();

	}
	void saveImageDCM(ITKImage3D::Pointer input, const char* path)
	{
		typedef itk::GDCMImageIO
			ImageIOType;
		ImageIOType::Pointer gdcmIO = ImageIOType::New();
		typedef  itk::ImageFileWriter<ITKImage3D> WriterType;
		WriterType::Pointer Writer = WriterType::New();
		std::stringstream ss;
		ss << path << ".dcm";
		Writer->SetImageIO(gdcmIO);
		Writer->SetFileName(ss.str().c_str());
		Writer->SetInput(input);
		Writer->Update();

	}
	//createCurveColony(int size) - It creates Colonies with radius from 1 to 'size'.
	//It puts in the global variable curveColony;
	void createCurveColony(int size){
		std::vector<std::vector<Crawler3D>> clusters;
		for (int i = 1; i <= size; i++){
			clusters = ConnectedElementsByRadius(i);
			curveColony.setElement(i, clusters.size());
		}

	}
	//createCurveScaleDistribuition(int size) - It counts the number of Clusters/Colony with radius from 1 to 'size'.
	//It puts in the global variable curveScaleDistribuition;
	void createCurveScaleDistribuition(int size){
		int cont = 0;
		for (int i = 1; i <= size; i++){
			for (int j = 0; j < global_clusters.size(); j++){
				if (global_clusters.at(j).size() == i) cont++;
			}
			curveScaleDistribuition.setElement(i, cont);
			cont = 0;
		}

	}

};
#endif