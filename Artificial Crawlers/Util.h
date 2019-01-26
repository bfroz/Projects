#ifndef Util_H
#define Util_H

#define PI 3.14146

#include <itkMinimumMaximumImageCalculator.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkGDCMImageIO.h>
#include <itkExtractImageFilter.h>
#include <tchar.h>
#include <windows.h>
#include <atlstr.h>
#include "ClassesUtil.h"
#include "getCirculosAneis.h"
#include <itkHoughTransform2DLinesImageFilter.h>
#include "itkThresholdImageFilter.h"
#include "itkRGBPixel.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkPNGImageIO.h"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

//Sismica
#include "itkDiscreteGaussianImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkScaleTransform.h"

typedef short PixelType;
typedef itk::Image<PixelType, 2> ITKImage2D;
typedef itk::Image<PixelType, 3> ITKImage3D;
typedef itk::ExtractImageFilter<ITKImage3D, ITKImage2D> ExtractImageFilter;
typedef  itk::ImageRegionIterator<ITKImage3D> Iterator3D;
typedef  itk::ImageRegionIterator<ITKImage2D> Iterator2D;
typedef itk::GDCMImageIO ImageIOType;
typedef itk::RGBPixel< unsigned char > RGBPixelType;
typedef itk::Image< RGBPixelType, 2 >    RGBImageType;

int countExamFolders(const char* pai) {

	TCHAR parent_path[200];
	_tcscpy_s(parent_path, CA2T(pai));
	// The hideous string manipulation code below
	// prepares a TCHAR wildcard string (sub_wild)
	// matching any subdirectory immediately under 
	// parent_path by appending "\*"
	size_t len = _tcslen(parent_path);
	const size_t alloc_len = len + 3;
	TCHAR* sub_wild = new TCHAR[alloc_len];
	_tcscpy_s(sub_wild, alloc_len, parent_path);
	if (len && sub_wild[len - 1] != _T('\\')) { sub_wild[len++] = _T('\\'); }
	sub_wild[len++] = _T('*');
	sub_wild[len++] = _T('\0');

	// File enumeration starts here
	WIN32_FIND_DATA fd;
	HANDLE hfind;
	int count = 0;
	if (INVALID_HANDLE_VALUE != (hfind = FindFirstFile(sub_wild, &fd))) {
		do {
			if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
				// is_alias_dir is true if directory name is "." or ".."
				const bool is_alias_dir = fd.cFileName[0] == _T('.') &&
					(!fd.cFileName[1] || (fd.cFileName[1] == _T('.') &&
					!fd.cFileName[2]));

				count += !is_alias_dir;
			}
		} while (FindNextFile(hfind, &fd));
		//CloseHandle(hfind);
	}

	delete[] sub_wild;
	return count;
}

ITKImage3D::Pointer openDicom(const char* path){
	typedef itk::ImageFileReader<ITKImage3D> ReaderType3D;
	ReaderType3D::Pointer itkReader = ReaderType3D::New();
	itkReader->SetFileName(path);
	ImageIOType::Pointer itkDicomIO = ImageIOType::New();
	itkReader->SetImageIO(itkDicomIO);
	itkReader->Update();
	return itkReader->GetOutput();
}

ITKImage3D::Pointer readImageDCM(const char* path){
	ITKImage3D::Pointer output;
	typedef itk::GDCMImageIO
		ImageIOType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	typedef  itk::ImageFileReader< ITKImage3D > ReaderType;
	try{
		ReaderType::Pointer reader = ReaderType::New();
		reader->SetImageIO(gdcmIO);
		std::stringstream ss;
		ss << path;
		reader->SetFileName(ss.str().c_str());
		reader->Update();
		output = reader->GetOutput();
	}
	catch (...){
		std::cout << "Couldn't read the image. Mind the path.\n";
	}
	return output;
}

void saveImageMHD(ITKImage2D::Pointer input, const char* path)
{
	typedef  itk::ImageFileWriter<ITKImage2D> WriterType;
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
void saveImageDCM(ITKImage2D::Pointer input, const char* path)
{
	typedef itk::GDCMImageIO
		ImageIOType;
	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	typedef  itk::ImageFileWriter<ITKImage2D> WriterType;
	WriterType::Pointer Writer = WriterType::New();
	std::stringstream ss;
	ss << path << ".dcm";
	Writer->SetImageIO(gdcmIO);
	Writer->SetFileName(ss.str().c_str());
	Writer->SetInput(input);
	Writer->Update();

}

float mean(std::vector<float> vetor){
	int cont = 0;
	for (int i = 0; i < vetor.size(); i++){
		cont += vetor.at(i);
	}
	return cont / vetor.size();
}

ITKImage3D::Pointer createImage3DWithBackground(int x, int y, int z, double background){
	ITKImage3D::Pointer output = ITKImage3D::New();
	ITKImage3D::SizeType size;
	size[0] = x;
	size[1] = y;
	size[2] = z;
	ITKImage3D::IndexType start;
	start.Fill(0);
	ITKImage3D::RegionType region(start, size);
	output->SetRegions(region);
	output->Allocate();
	output->FillBuffer(background);

	return output;
}

double radian_to_angle(double radian){

	return (radian * 180) / PI;
}

double mean_direction(std::vector<float> vetores){

	double sum_seno = 0;
	double sum_cos = 0;

	double rho;

	//radians = angle * PI/180;

	double radian;
	for (std::vector<float>::iterator it = vetores.begin(); it != vetores.end(); ++it){
		radian = *it;

		sum_seno = sum_seno + sin(radian);
		sum_cos = sum_cos + cos(radian);


	}
	if (sum_seno > 0 && sum_cos > 0) rho = 0;
	else rho = PI;

	return atan(sum_seno / sum_cos + rho);
}

double sample_circular_variance(double mean_direction){

	return 1 - mean_direction;

}

double sample_circular_standard_deviation(double mean_direction){
	if (mean_direction < 0) mean_direction = mean_direction *-1;

	double result = -2 * log(mean_direction);
	if (result < 0) result = result *-1;
	return sqrt(result);
}

double resultant_vector_strenght(std::vector<float> vetores){
	double sum_seno = 0;
	double sum_cos = 0;

	//radians = angle * PI/180;

	double radian;
	for (std::vector<float>::iterator it = vetores.begin(); it != vetores.end(); ++it){
		radian = *it;

		sum_seno = sum_seno + sin(radian);
		sum_cos = sum_cos + cos(radian);
	}

	return sqrt(pow(sum_seno, 2) + pow(sum_cos, 2)) / vetores.size();
}

double kurtosis(std::vector<float> vetores, double mean_direction){
	double sum_cos = 0;
	double radian;
	for (std::vector<float>::iterator it = vetores.begin(); it != vetores.end(); ++it){
		radian = *it;

		sum_cos = sum_cos + cos(2 * (radian - mean_direction));
	}

	return (1.0 / vetores.size()) * sum_cos;

}

double skewness(std::vector<float> vetores, double mean_direction){
	double sum_sin = 0;
	double radian;
	for (std::vector<float>::iterator it = vetores.begin(); it != vetores.end(); ++it){
		radian = *it;

		sum_sin = sum_sin + sin(2 * (*it - mean_direction));
	}
	return (1.0 / vetores.size()) * sum_sin;

}

double angle_to_radian(double angle){

	return (angle * PI) / 180;
}

ITKImage2D::Pointer extrairFatia(ITKImage3D::Pointer image, int index){
	ExtractImageFilter::Pointer slicingFilter = ExtractImageFilter::New();

	ITKImage3D::RegionType inputRegion = image->GetLargestPossibleRegion();
	ITKImage3D::SizeType size = inputRegion.GetSize();
	size[2] = 0;

	ITKImage3D::IndexType start = inputRegion.GetIndex();
	start[2] = index;

	ITKImage3D::RegionType desiredRegion;
	desiredRegion.SetSize(size);
	desiredRegion.SetIndex(start);
	slicingFilter->SetExtractionRegion(desiredRegion);
	slicingFilter->SetInput(image);

	slicingFilter->Update();

	return slicingFilter->GetOutput();
}

//sobelGradientCalculator(...) - Calculates the Sobel Gradient of the neighborhood elements
float sobelGradientCalculator(NeighborhoodIterator it, int raio){
	NeighborhoodInnerProduct IP;
	itk::Size<2> radius;
	radius.Fill(raio);
	itk::SobelOperator<PixelType, 2> sH;
	sH.SetDirection(0);
	sH.CreateToRadius(radius);
	sH.CreateDirectional();
	itk::SobelOperator<PixelType, 2> sV;
	sV.SetDirection(1);
	sV.CreateToRadius(radius);
	sV.CreateDirectional();
	double Gx, Gy;
	Gx = IP(it, sH);
	Gy = IP(it, sV);
	float direction = atan2(Gy, Gx);
	return direction;
}

float sobelGradientX(NeighborhoodIterator it, int raio){
	NeighborhoodInnerProduct IP;
	itk::Size<2> radius;
	radius.Fill(raio);
	itk::SobelOperator<PixelType, 2> sH;
	sH.SetDirection(0);
	sH.CreateToRadius(radius);
	sH.CreateDirectional();
	double Gx;
	Gx = IP(it, sH);
	return Gx;
}

float sobelGradientY(NeighborhoodIterator it, int raio){
	NeighborhoodInnerProduct IP;
	itk::Size<2> radius;
	radius.Fill(raio);
	itk::SobelOperator<PixelType, 2> sV;
	sV.SetDirection(1);
	sV.CreateToRadius(radius);
	sV.CreateDirectional();
	double Gy;
	Gy = IP(it, sV);
	return Gy;
}

float totalNonBackgroundPixels(ITKImage2D::Pointer image){
	typedef itk::MinimumMaximumImageCalculator <ITKImage2D> ImageCalculatorFilterType;
	typedef  itk::ImageRegionIterator<ITKImage2D> Iterator;

	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->ComputeMinimum();
	float minimum = imageCalculatorFilter->GetMinimum();

	Iterator  it(image, image->GetLargestPossibleRegion());
	int cont = 0;
	for (it = it.Begin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() > minimum) cont++;
	}
	return cont;
}

int getMinimum2D(ITKImage2D::Pointer image){
	typedef itk::MinimumMaximumImageCalculator <ITKImage2D> ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->ComputeMinimum();
	return imageCalculatorFilter->GetMinimum();

}
int getMinimum3D(ITKImage3D::Pointer image){
	typedef itk::MinimumMaximumImageCalculator <ITKImage3D> ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->ComputeMinimum();
	return imageCalculatorFilter->GetMinimum();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											//Funções de treino dos arquivos das curvas
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Curve_T GerarTemplateMedio(std::vector<Curve_T> candidatos, int dominio){
	Curve_T temp;
	int size = candidatos.size();
	//Media das curvas (Regressão)
	float x = 0, y = 0;
	for (int i = 0; i < dominio; i++){
		for (int j = 0; j < size; j++){
			x += candidatos.at(j).at(i).x;
			y += candidatos.at(j).at(i).y;
		}
		temp.push_back(PointN(x / size, y / size));
		x = 0;
		y = 0;
	}
	return temp;
}

Curve_T InterpretarCurva(const char* diretorio){
	std::string line;
	std::string token;
	std::ifstream file(diretorio);
	std::vector<float> coords;
	Curve_T curve;
	if (file.is_open())
	{
		while (getline(file, line))
		{
			std::istringstream ss(line);
			while (std::getline(ss, token, ',')){
				coords.push_back(atof(token.c_str()));
			}
			if (coords.size() > 0) curve.push_back(PointN(coords[0], coords[1]));
			coords.clear();
		}
		file.close();
	}
	return curve;
}

float InterpretarArea(const char* diretorio){
	std::string line;
	std::string token;
	std::ifstream file(diretorio);
	float area;
	Curve_T curve;
	if (file.is_open())
	{
		while (getline(file, line))
		{
			area = atof(line.c_str());
			break;
		}
		file.close();
	}
	return area;
}
std::vector<ITKImage3D::Pointer> encontrarAneisEsferas(ITKImage3D::Pointer img3D){
	VOI* voi = new VOI();
	LungSegmentedRegion * lungSegReg = new LungSegmentedRegion(img3D, voi, new ITKInformationList());
	getCirculosAneis * esferas_aneis = new getCirculosAneis(lungSegReg);
	vector<ITKImage3D::Pointer>* esfera_aneis_vector;
	vector<ITKImage3D::Pointer> aneis_esferas_vector;
	esfera_aneis_vector = esferas_aneis->executa();

	for (int i = 1; i < esfera_aneis_vector->size(); i++){
		aneis_esferas_vector.push_back(esfera_aneis_vector->at(i));
	}
	return aneis_esferas_vector;
}

void salvarAneisEsferas(std::vector<ITKImage3D::Pointer> elementos, const char * diretorio){
	typedef  itk::ImageFileWriter< ITKImage3D  > WriterType;
	WriterType::Pointer writer = WriterType::New();
	std::stringstream ssEsfera;
	std::stringstream ssAnel;
	for (int i = 0; i < 5; i++){
		ssEsfera << diretorio << "esfera" << i + 1 << ".dcm";
		writer->SetFileName(ssEsfera.str());
		writer->SetInput(elementos.at(i));
		writer->Update();
		ssEsfera.str("");
	}
	for (int i = 6; i < elementos.size(); i++){
		ssAnel << diretorio << "anel" << i - 5 << ".dcm";
		writer->SetFileName(ssAnel.str());
		writer->SetInput(elementos.at(i));
		writer->Update();
		ssAnel.str("");
	}
}


////SISMICA
//3D Functions
ITKImage3D::Pointer GaussianITK(ITKImage3D::Pointer input, float variance){
	typedef itk::DiscreteGaussianImageFilter<ITKImage3D, ITKImage3D>  filterType;
	filterType::Pointer gaussianFilter = filterType::New();
	gaussianFilter->SetInput(input);
	gaussianFilter->SetVariance(variance);
	gaussianFilter->Update();
	return gaussianFilter->GetOutput();
}

ITKImage3D::Pointer EQHistograma(ITKImage3D::Pointer input){

	int x = input->GetLargestPossibleRegion().GetSize()[0];

	int y = input->GetLargestPossibleRegion().GetSize()[1];

	int z = input->GetLargestPossibleRegion().GetSize()[2];


	typedef itk::MinimumMaximumImageCalculator< ITKImage3D > MinimumMaximumImageCalculatorType;

	MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter = MinimumMaximumImageCalculatorType::New();
	minimumMaximumImageCalculatorFilter->SetImage(input);
	minimumMaximumImageCalculatorFilter->Compute();

	short min = minimumMaximumImageCalculatorFilter->GetMinimum();
	short max = minimumMaximumImageCalculatorFilter->GetMaximum();

	short size;
	if (min < 0) size = -2 * min;
	else size = 2 * min;

	std::vector<int> histo;
	//unsigned int* histo = new unsigned int[size];



	for (int i = 0; i < size; i++) {
		histo[i] = 0;
	}

	for (int i = 0; i < x; i++){
		for (int j = 0; j < y; j++){
			for (int k = 0; k < z; k++){

				itk::Index<3> center;
				long centerId[] = { i, j };
				center.SetIndex(centerId);

				short val = input->GetPixel(center);

				if (val != min){
					int id = (val)-(min);


					histo[id]++;
				}
			}
		}
	}

	for (int i = 1; i < size; i++){
		histo[i] = histo[i - 1] + histo[i];
	}

	unsigned int* histoNovo = new unsigned int[size];

	for (int i = 0; i < size; i++){
		double numerador = size;
		double denominador = x*y;
		double histoVal = histo[i];
		double val = (numerador / denominador)*histoVal;

		histoNovo[i] = val;
	}

	for (int i = 0; i < x; i++){
		for (int j = 0; j < y; j++){
			for (int k = 0; k < z; k++){

				itk::Index<3> center;
				long centerId[] = { i, j };
				center.SetIndex(centerId);

				short val = input->GetPixel(center);

				short id = (val)-(min);

				if (id == 0)
					val = min;
				else
					val = histoNovo[id];

				input->SetPixel(center, val);



			}
		}
	}
	//delete histo;
	delete histoNovo;
	return input;
}

//OtsuITK(...) - Runs Otsu Algorithm. 
//- outside - Value of the output higher threshold
//- outside - Value of the output lower threshold
//- hist - Number of Histogram Bins
ITKImage3D::Pointer OtsuITK(ITKImage3D::Pointer input, int outside, int inside, int hist){
	typedef itk::OtsuThresholdImageFilter<ITKImage3D, ITKImage3D>  FilterType;
	FilterType::Pointer OtsuFilter = FilterType::New();
	OtsuFilter->SetInput(input);
	OtsuFilter->SetOutsideValue(outside);
	OtsuFilter->SetInsideValue(inside);
	OtsuFilter->SetNumberOfHistogramBins(hist);
	OtsuFilter->Update();
	return OtsuFilter->GetOutput();
}

//invertImage(...) - Inverts a 3D binary ITKImage3D
ITKImage3D::Pointer invertImage(ITKImage3D::Pointer input){
	ITKImage3D::Pointer output = ITKImage3D::New();
	ITKImage3D::SizeType size;
	size[0] = input->GetLargestPossibleRegion().GetSize()[0];
	size[1] = input->GetLargestPossibleRegion().GetSize()[1];
	size[2] = input->GetLargestPossibleRegion().GetSize()[2];
	ITKImage3D::IndexType start;
	start.Fill(0);
	ITKImage3D::RegionType region(start, size);
	output->SetRegions(region);
	output->Allocate();

	Iterator3D  it(input, input->GetLargestPossibleRegion());
	for (it = it.Begin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() > 0) output->SetPixel(it.GetIndex(), 0);
		else output->SetPixel(it.GetIndex(), 255);
	}
	return output;
}


//2D Functions
ITKImage2D::Pointer GaussianITK(ITKImage2D::Pointer input, float variance){
	typedef itk::DiscreteGaussianImageFilter<ITKImage2D, ITKImage2D>  filterType;
	filterType::Pointer gaussianFilter = filterType::New();
	gaussianFilter->SetInput(input);
	gaussianFilter->SetVariance(variance);
	gaussianFilter->Update();
	return gaussianFilter->GetOutput();
}

ITKImage2D::Pointer EQHistograma(ITKImage2D::Pointer input){

	int x = input->GetLargestPossibleRegion().GetSize()[0];

	int y = input->GetLargestPossibleRegion().GetSize()[1];

	//int z = input->GetLargestPossibleRegion().GetSize()[2];
	int z = 1;


	typedef itk::MinimumMaximumImageCalculator< ITKImage2D > MinimumMaximumImageCalculatorType;

	MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter = MinimumMaximumImageCalculatorType::New();
	minimumMaximumImageCalculatorFilter->SetImage(input);
	minimumMaximumImageCalculatorFilter->Compute();

	short min = minimumMaximumImageCalculatorFilter->GetMinimum();
	short max = minimumMaximumImageCalculatorFilter->GetMaximum();

	short size;
	if (min < 0) size = -2 * min;
	else size = 2 * min;

	//std::vector<int> histo;
	unsigned int* histo = new unsigned int[size];



	for (int i = 0; i < size; i++) {
		histo[i] = 0;
	}

	for (int i = 0; i < x; i++){
		for (int j = 0; j < y; j++){
			for (int k = 0; k < z; k++){

				itk::Index<2> center;
				long centerId[] = { i, j };
				center.SetIndex(centerId);

				short val = input->GetPixel(center);

				if (val != min){
					int id = (val)-(min);


					histo[id]++;
				}
			}
		}
	}

	for (int i = 1; i < size; i++){
		histo[i] = histo[i - 1] + histo[i];
	}

	unsigned int* histoNovo = new unsigned int[size];

	for (int i = 0; i < size; i++){
		double numerador = size;
		double denominador = x*y;
		double histoVal = histo[i];
		double val = (numerador / denominador)*histoVal;

		histoNovo[i] = val;
	}

	for (int i = 0; i < x; i++){
		for (int j = 0; j < y; j++){
			for (int k = 0; k < z; k++){

				itk::Index<2> center;
				long centerId[] = { i, j };
				center.SetIndex(centerId);

				short val = input->GetPixel(center);

				short id = (val)-(min);

				if (id == 0)
					val = min;
				else
					val = histoNovo[id];

				input->SetPixel(center, val);



			}
		}
	}
	//delete histo;
	delete histoNovo;
	return input;
}

ITKImage2D::Pointer HistogramEqualizationITK(ITKImage2D::Pointer input, int radius){
	typedef  itk::AdaptiveHistogramEqualizationImageFilter< ITKImage2D > AdaptiveHistogramEqualizationImageFilterType;
	AdaptiveHistogramEqualizationImageFilterType::Pointer adaptiveHistogramEqualizationImageFilter = AdaptiveHistogramEqualizationImageFilterType::New();
	ITKImage2D::SizeType size;
	size.Fill(radius);
	adaptiveHistogramEqualizationImageFilter->SetInput(input);
	adaptiveHistogramEqualizationImageFilter->SetRadius(size);
	adaptiveHistogramEqualizationImageFilter->Update();
	return adaptiveHistogramEqualizationImageFilter->GetOutput();
}

//OtsuITK(...) - Runs Otsu Algorithm. 
//- outside - Value of the output higher threshold
//- outside - Value of the output lower threshold
//- hist - Number of Histogram Bins
ITKImage2D::Pointer OtsuITK(ITKImage2D::Pointer input, int outside, int inside, int hist){
	typedef itk::OtsuThresholdImageFilter<ITKImage2D, ITKImage2D>  FilterType;
	FilterType::Pointer OtsuFilter = FilterType::New();
	OtsuFilter->SetInput(input);
	OtsuFilter->SetOutsideValue(outside);
	OtsuFilter->SetInsideValue(inside);
	OtsuFilter->SetNumberOfHistogramBins(hist);
	OtsuFilter->Update();
	return OtsuFilter->GetOutput();
}

//invertImage(...) - Inverts a 3D binary ITKImage3D
ITKImage2D::Pointer invertImage(ITKImage2D::Pointer input){
	ITKImage2D::Pointer output = ITKImage2D::New();
	ITKImage2D::SizeType size;
	size[0] = input->GetLargestPossibleRegion().GetSize()[0];
	size[1] = input->GetLargestPossibleRegion().GetSize()[1];
	ITKImage2D::IndexType start;
	start.Fill(0);
	ITKImage2D::RegionType region(start, size);
	output->SetRegions(region);
	output->Allocate();

	Iterator2D  it(input, input->GetLargestPossibleRegion());
	for (it = it.Begin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() > 0) output->SetPixel(it.GetIndex(), 0);
		else output->SetPixel(it.GetIndex(), 255);
	}
	return output;
}
ITKImage2D::Pointer Erode(ITKImage2D::Pointer input, int kernelx, int kernely){
	ITKImage2D::Pointer image = input;
	ITKImage2D::SizeType size;
	size[0] = input->GetLargestPossibleRegion().GetSize()[0];
	size[1] = input->GetLargestPossibleRegion().GetSize()[1];
	for (int i = 0; i< kernelx; i++){
		ITKImage2D::Pointer output = ITKImage2D::New();
		const unsigned int numberOfPixels = size[0] * size[1];
		ITKImage2D::IndexType start;
		start.Fill(0);
		ITKImage2D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		itk::ImageRegionIterator<ITKImage2D> out(output, region);
		NeighborhoodIterator::RadiusType raio;
		for (unsigned int i = 0; i < ITKImage2D::ImageDimension; ++i) raio[i] = 1;
		NeighborhoodIterator it(raio, image, image->GetLargestPossibleRegion());
		int central = it.Size() / 2;
		int right_pixel, left_pixel, below_pixel, above_pixel;
		int s = it.GetStride(1);

		//Na direção X
			out = out.Begin();
		it.GoToBegin();
		while (!it.IsAtEnd()){
			right_pixel = it.GetPixel(central + 1);
			left_pixel = it.GetPixel(central - 1);
			if (right_pixel > 0 && left_pixel > 0 && it.GetCenterPixel() > 0){
				out.Set(255);
			}
			else{
				out.Set(0);
			}
			++it;
			++out;
		}
		image = output;
	}
	for (int i = 0; i< kernely; i++){
		ITKImage2D::Pointer output = ITKImage2D::New();
		const unsigned int numberOfPixels = size[0] * size[1];
		ITKImage2D::IndexType start;
		start.Fill(0);
		ITKImage2D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		itk::ImageRegionIterator<ITKImage2D> out(output, region);
		NeighborhoodIterator::RadiusType raio;
		for (unsigned int i = 0; i < ITKImage2D::ImageDimension; ++i) raio[i] = 1;
		NeighborhoodIterator it(raio, image, image->GetLargestPossibleRegion());
		int central = it.Size() / 2;
		int right_pixel, left_pixel, below_pixel, above_pixel;
		int s = it.GetStride(1);
		//Na direção de Y
			out = out.Begin();
		it.GoToBegin();
		while (!it.IsAtEnd()){
			above_pixel = it.GetPixel(central - s);
			below_pixel = it.GetPixel(central + s);
			if (above_pixel > 0 && below_pixel > 0 && it.GetCenterPixel() > 0){
				out.Set(255);
			}
			else{
				out.Set(0);
			}
			++it;
			++out;
		}
		image = output;
	}

	return image;
}
ITKImage2D::Pointer Dilate(ITKImage2D::Pointer input, int kernelx, int kernely){
	ITKImage2D::Pointer image = input;
	ITKImage2D::SizeType size;
	size[0] = input->GetLargestPossibleRegion().GetSize()[0];
	size[1] = input->GetLargestPossibleRegion().GetSize()[1];

	for (int i = 0; i< kernelx; i++){
		ITKImage2D::Pointer output = ITKImage2D::New();
		const unsigned int numberOfPixels = size[0] * size[1];
		ITKImage2D::IndexType start;
		start.Fill(0);
		ITKImage2D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		itk::ImageRegionIterator<ITKImage2D> out(output, region);
		NeighborhoodIterator::RadiusType raio;
		for (unsigned int i = 0; i < ITKImage2D::ImageDimension; ++i) raio[i] = 1;
		NeighborhoodIterator it(raio, image, image->GetLargestPossibleRegion());
		int central = it.Size() / 2;
		int right_pixel, left_pixel, below_pixel, above_pixel;
		int s = it.GetStride(1);
		//Na direção X
			out = out.Begin();
		it.GoToBegin();
		while (!it.IsAtEnd()){
			right_pixel = it.GetPixel(central + 1);
			left_pixel = it.GetPixel(central - 1);
			if (left_pixel == 0 && right_pixel == 0 && it.GetCenterPixel() == 0){
				out.Set(0);
			}
			else{
				out.Set(255);
			}
			++it;
			++out;
		}
		image = output;
	}
	for (int i = 0; i< kernely; i++){
		ITKImage2D::Pointer output = ITKImage2D::New();
		const unsigned int numberOfPixels = size[0] * size[1];
		ITKImage2D::IndexType start;
		start.Fill(0);
		ITKImage2D::RegionType region(start, size);
		output->SetRegions(region);
		output->Allocate();
		itk::ImageRegionIterator<ITKImage2D> out(output, region);
		NeighborhoodIterator::RadiusType raio;
		for (unsigned int i = 0; i < ITKImage2D::ImageDimension; ++i) raio[i] = 1;
		NeighborhoodIterator it(raio, image, image->GetLargestPossibleRegion());
		int central = it.Size() / 2;
		int right_pixel, left_pixel, below_pixel, above_pixel;
		int s = it.GetStride(1);
		//Na direção de Y
			out = out.Begin();
		it.GoToBegin();
		while (!it.IsAtEnd()){
			above_pixel = it.GetPixel(central - s);
			below_pixel = it.GetPixel(central + s);
			if (above_pixel == 0 && below_pixel == 0 && it.GetCenterPixel() == 0){
				out.Set(0);
			}
			else{
				out.Set(255);
			}
			++it;
			++out;
		}
		image = output;
	}

	return image;
}

ITKImage2D::Pointer or(ITKImage2D::Pointer a, ITKImage2D::Pointer b){
	ITKImage2D::Pointer output = ITKImage2D::New();
	ITKImage2D::SizeType size;
	size[0] = a->GetLargestPossibleRegion().GetSize()[0];
	size[1] = a->GetLargestPossibleRegion().GetSize()[1];
	ITKImage2D::IndexType start;
	start.Fill(0);
	ITKImage2D::RegionType region(start, size);
	output->SetRegions(region);
	output->Allocate();

	Iterator2D  it(a, a->GetLargestPossibleRegion());
	for (it = it.Begin(); !it.IsAtEnd(); ++it)
	{
		if (it.Get() > 0 || b->GetPixel(it.GetIndex()) > 0)
			output->SetPixel(it.GetIndex(), 255);
	}
	
	return output;
}

ITKImage2D::Pointer houghLineITK(ITKImage2D::Pointer image, int num_lines, float radius, float variance){
	typedef itk::ThresholdImageFilter<ITKImage2D> ThresholdFilterType;
	ThresholdFilterType::Pointer threshFilter = ThresholdFilterType::New();
	threshFilter->SetInput(image);
	threshFilter->SetOutsideValue(0);
	unsigned char threshBelow = 0;
	unsigned char threshAbove = 255;
	threshFilter->ThresholdOutside(threshBelow, threshAbove);
	threshFilter->Update();

	typedef  itk::HoughTransform2DLinesImageFilter< PixelType, PixelType > HoughTransform2DLinesImageFilter;
	HoughTransform2DLinesImageFilter::Pointer houghTransform2DLinesImageFilter = HoughTransform2DLinesImageFilter::New();
	houghTransform2DLinesImageFilter->SetInput(threshFilter->GetOutput());
	houghTransform2DLinesImageFilter->SetNumberOfLines(num_lines);
	houghTransform2DLinesImageFilter->SetDiscRadius(radius);
	houghTransform2DLinesImageFilter->SetVariance(variance);
	houghTransform2DLinesImageFilter->Update();

	ITKImage2D::Pointer localAccumulator = houghTransform2DLinesImageFilter->GetOutput();
	HoughTransform2DLinesImageFilter::LinesListType lines;
	lines = houghTransform2DLinesImageFilter->GetLines(num_lines);

	typedef  unsigned char                            OutputPixelType;
	typedef  itk::Image< OutputPixelType, 2 > OutputITKImage2D;
	OutputITKImage2D::Pointer  localOutputImage = OutputITKImage2D::New();
	OutputITKImage2D::RegionType region;
	region.SetSize(image->GetLargestPossibleRegion().GetSize());
	region.SetIndex(image->GetLargestPossibleRegion().GetIndex());
	localOutputImage->SetRegions(region);
	localOutputImage->SetOrigin(image->GetOrigin());
	localOutputImage->SetSpacing(image->GetSpacing());
	localOutputImage->Allocate();

	typedef HoughTransform2DLinesImageFilter::LinesListType::const_iterator LineIterator;
	LineIterator itLines = lines.begin();
	while (itLines != lines.end())
	{
		typedef HoughTransform2DLinesImageFilter::LineType::PointListType  PointListType;
		PointListType pointsList = (*itLines)->GetPoints();
		PointListType::const_iterator   itPoints = pointsList.begin();
		double u[2];
		u[0] = (*itPoints).GetPosition()[0];
		u[1] = (*itPoints).GetPosition()[1];
		itPoints++;
		double v[2];
		v[0] = u[0] - (*itPoints).GetPosition()[0];
		v[1] = u[1] - (*itPoints).GetPosition()[1];
		double norm = std::sqrt(v[0] * v[0] + v[1] * v[1]);
		v[0] /= norm;
		v[1] /= norm;
		OutputITKImage2D::IndexType localIndex;
		itk::Size<2> size = localOutputImage->GetLargestPossibleRegion().GetSize();
		float diag = std::sqrt((float)(size[0] * size[0] + size[1] * size[1]));
		for (int i = static_cast<int>(-diag); i<static_cast<int>(diag); i++)
		{
			localIndex[0] = (long int)(u[0] + i*v[0]);
			localIndex[1] = (long int)(u[1] + i*v[1]);
			OutputITKImage2D::RegionType outputRegion =
				localOutputImage->GetLargestPossibleRegion();
			if (outputRegion.IsInside(localIndex))
			{
				localOutputImage->SetPixel(localIndex, 255);
			}
		}
		itLines++;
	}
	typedef itk::CastImageFilter< OutputITKImage2D, ITKImage2D > CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(localOutputImage);
	castFilter->Update();

	return castFilter->GetOutput();
}

ITKImage2D::Pointer RescaleITK(ITKImage2D::Pointer input, int minimum, int maximum){
	typedef itk::RescaleIntensityImageFilter< ITKImage2D, ITKImage2D > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(input);
	rescaleFilter->SetOutputMinimum(minimum);
	rescaleFilter->SetOutputMaximum(maximum);
	rescaleFilter->Update();
	return rescaleFilter->GetOutput();
}
RGBImageType::Pointer paintFaults(ITKImage2D::Pointer input, ITKImage2D::Pointer binary){
	int sizex = input->GetLargestPossibleRegion().GetSize()[0];
	int sizey = input->GetLargestPossibleRegion().GetSize()[1];

	typedef itk::CastImageFilter< ITKImage2D, RGBImageType > CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(RescaleITK(input, 0, 255));
	castFilter->Update();

	RGBImageType::Pointer output = castFilter->GetOutput();
	RGBImageType::PixelType RedPixel;
	RedPixel.SetRed(255);
	RedPixel.SetGreen(0);
	RedPixel.SetBlue(0);

	RGBImageType::IndexType ColorIndex;
	ITKImage2D::IndexType index;
	for (int i = 0; i< sizex; i++){
		for (int j = 0; j < sizey; j++){
			index[0] = i;
			index[1] = j;
			ColorIndex[0] = i;
			ColorIndex[1] = j;
			if (binary->GetPixel(index) == 255){
				output->SetPixel(ColorIndex, RedPixel);
			}
		}
	}
	return output;
}

void saveImagePNG(RGBImageType::Pointer input, const char* path)
{
	typedef itk::PNGImageIO ImageIOType;
	ImageIOType::Pointer pngIO = ImageIOType::New();

	/*typedef itk::CastImageFilter< RGBImageType, unsigned char > CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(input);
	castFilter->Update();*/

	typedef  itk::ImageFileWriter<RGBImageType> WriterType;
	WriterType::Pointer Writer = WriterType::New();
	std::stringstream ss;
	ss << path << ".png";
	Writer->SetImageIO(pngIO);
	Writer->SetFileName(ss.str().c_str());
	Writer->SetInput(input);
	Writer->Update();

}
void saveFaults(std::vector<RGBImageType::Pointer> faults, const char * diretorio){
	std::stringstream ss;
	
	for (int i = 0; i < faults.size(); i++){
		ss << diretorio << "\\" << i;
		saveImagePNG(faults.at(i), ss.str().c_str());
		ss.str("");
	}
}
//void RedimensionarPNG(const char * path){
//	typedef itk::PNGImageIO ImageIOType;
//	ImageIOType::Pointer pngIO = ImageIOType::New();
//	typedef itk::ImageFileReader<RGBImageType> ReaderType2D;
//	ReaderType2D::Pointer itkReader = ReaderType2D::New();
//	itkReader->SetFileName(path);
//	itkReader->SetImageIO(pngIO);
//	itkReader->Update();
//	//TODO RESAMPLE
//	typedef double                              ScalarType;
//	RGBImageType::Pointer inputImage = itkReader->GetOutput();
//
//	RGBImageType::RegionType region = inputImage->GetLargestPossibleRegion();
//	RGBImageType::SizeType size = region.GetSize();
//	RGBImageType::SpacingType spacing = inputImage->GetSpacing();
//
//	itk::Index< 2 > centralPixel;
//	centralPixel[0] = size[0] / 2;
//	centralPixel[1] = size[1] / 2;
//	itk::Point< ScalarType, 2 > centralPoint;
//	centralPoint[0] = centralPixel[0];
//	centralPoint[1] = centralPixel[1];
//
//	typedef itk::ScaleTransform< ScalarType, 2 > ScaleTransformType;
//	ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();
//
//	ScaleTransformType::ParametersType parameters = scaleTransform->GetParameters();
//	parameters[0] = 3;
//	parameters[1] = 3;
//
//	scaleTransform->SetParameters(parameters);
//	scaleTransform->SetCenter(centralPoint);
//
//	typedef itk::LinearInterpolateImageFunction< RGBImageType, ScalarType > LinearInterpolatorType;
//	LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
//
//	typedef itk::ResampleImageFilter< RGBImageType, RGBImageType > ResampleFilterType;
//	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
//
//	resampleFilter->SetInput(inputImage);
//	resampleFilter->SetTransform(scaleTransform);
//	resampleFilter->SetInterpolator(interpolator);
//	resampleFilter->SetSize(size);
//	resampleFilter->SetOutputSpacing(spacing);
//
//	std::stringstream ss;
//	ss << path << "-r.png";
//
//	saveImagePNG(resampleFilter->GetOutput(), ss.str().c_str());
//}

ITKImage2D::Pointer SobelImage(ITKImage2D::Pointer image, int direction){
	typedef itk::MinimumMaximumImageCalculator <ITKImage2D> ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->ComputeMinimum();
	float minimum = imageCalculatorFilter->GetMinimum();

	ITKImage2D::Pointer  localOutputImage = ITKImage2D::New();
	ITKImage2D::RegionType region;
	region.SetSize(image->GetLargestPossibleRegion().GetSize());
	region.SetIndex(image->GetLargestPossibleRegion().GetIndex());
	localOutputImage->SetRegions(region);
	/*localOutputImage->SetOrigin(image->GetOrigin());
	localOutputImage->SetSpacing(image->GetSpacing());*/
	localOutputImage->Allocate();
	localOutputImage->FillBuffer(minimum);
	NeighborhoodIterator::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	float value;
	NeighborhoodIterator it(radius, image, image->GetLargestPossibleRegion());
	it.GoToBegin();
	/*NeighborhoodIterator it2(radius, localOutputImage, localOutputImage->GetLargestPossibleRegion());
	it2.GoToBegin();*/
	while (!it.IsAtEnd())
	{
		if (direction == 0){
			value = sobelGradientX(it, 1);
		}
		else{
			value = sobelGradientY(it, 1);
		}
		if (value == 0) localOutputImage->SetPixel(it.GetIndex(), minimum);
		else localOutputImage->SetPixel(it.GetIndex(), value);
		

		++it;
		/*++it2;*/
	}
	return localOutputImage;
}
#endif