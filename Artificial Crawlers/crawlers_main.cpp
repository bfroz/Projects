#include "Util.h"
#include "ArtificialCrawler3D.h"
#include <cstdlib>

int main(){
	std::vector<double> caracteristicas; //Vetor de caracteristicas

	ITKImage3D::Pointer imagem = readImageDCM("~PATH~/imagem.dcm");
	ArtificialCrawler3D acrawlers = ArtificialCrawler3D(imagem, 5); //Coloco o valor 5 como valor de energia inicial de cada crawler;
	acrawlers.start(); //Essa função de start() sem parâmetro roda quantas vezes for necessário os ciclos do ACrawlers, com critério de parada sendo quando nenhum indivíduo se moveu;

	Curve_T HabitantSettlement = acrawlers.getCurveHabitantsSettlement(); //Esse tipo "Curve_T" é um vector<PointN>
	Curve_T ColonyFormation = acrawlers.getCurveColonyFormation();
	Curve_T NumberCrawlers = acrawlers.getCurveNumberCrawlers();
	Curve_T ScaleDistribution = acrawlers.getCurveScaleDistribution();
	//Essas variáveis com "tem" no começo são os templates. Exemplo, temColonyFormation ou temNumberCrawlers.
	//Elas também são elementos do tipo Curve_T.
	//Abstraí como conseguir eles, por que é igual o feito acima. Se quiser criar um template, de preferência, 
	//chame um procedimento antes que faça o mesmo que essa main aqui faz e passe como parâmetro (ou global) as curvas depois.
	
	//Euclidian
	caracteristicas.push_back(curve_distance(ColonyFormation, temColonyFormation)); //Colony Formation
	caracteristicas.push_back(curve_distance(NumberCrawlers, temNumberCrawlers)); // Number of Crawlers
	caracteristicas.push_back(curve_distance(ScaleDistribution, temScaleDistribution)); // Scale Distribution
	caracteristicas.push_back(curve_distance(HabitantSettlement, temHabitantSettlement)); // Habitants Settlement

	//Jaccard
	caracteristicas.push_back(curve_jaccard_distance(ColonyFormation, temColonyFormation, 0.1)); //Colony Formation
	caracteristicas.push_back(curve_jaccard_distance(NumberCrawlers, temNumberCrawlers, 0.1)); // Number of Crawlers
	caracteristicas.push_back(curve_jaccard_distance(ScaleDistribution, temScaleDistribution, 0.1)); // Scale Distribution
	caracteristicas.push_back(curve_jaccard_distance(HabitantSettlement, temHabitantSettlement, 0.1)); // Habitants Settlement

	//Simple Matching
	caracteristicas.push_back(curve_simple_matching(ColonyFormation, temColonyFormation, 0.1)); //Colony Formation
	caracteristicas.push_back(curve_simple_matching(NumberCrawlers, temNumberCrawlers, 0.1)); // Number of Crawlers
	caracteristicas.push_back(curve_simple_matching(ScaleDistribution, temScaleDistribution, 0.1)); // Scale Distribution
	caracteristicas.push_back(curve_simple_matching(HabitantSettlement, temHabitantSettlement, 0.1)); // Habitants Settlement

	//Manhattan
	caracteristicas.push_back(curve_manhattan_distance(ColonyFormation, temColonyFormation, 0.1)); //Colony Formation
	caracteristicas.push_back(curve_manhattan_distance(NumberCrawlers, temNumberCrawlers, 0.1)); // Number of Crawlers
	caracteristicas.push_back(curve_manhattan_distance(ScaleDistribution, temScaleDistribution, 0.1)); // Scale Distribution
	caracteristicas.push_back(curve_manhattan_distance(HabitantSettlement, temHabitantSettlement, 0.1)); // Habitants Settlement

	//Areas

	double curve_colony = acrawlers.getCurveColonyFormationnArea();
	double curve_number = acrawlers.getCurveNumberCrawlersArea();
	double curve_distribution = acrawlers.getCurveScaleDistributionArea();
	double curve_habitants = acrawlers.getCurveHabitantsSettlementArea();
	caracteristicas.push_back(curve_colony); //Colony Formation
	caracteristicas.push_back(curve_number); // Number of Crawlers
	caracteristicas.push_back(curve_distribution); // Scale Distribution
	caracteristicas.push_back(curve_habitants); // Habitants Settlement

	return 0;
}
