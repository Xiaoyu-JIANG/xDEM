#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>
#include <iostream>

#include "src/QuadronAnalysis.h"

using namespace std;

#pragma warning(disable : 4996)

void oneCase(const string path)
{
	cout << path << endl;

	// Create output files for quadron analysis
	std::stringstream fileName;
	fileName << path << "quadronAnalysis_cell.dat";
	FILE* fileCell = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_particle.dat";
	FILE* fileParticle = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_quadron.dat";
	FILE* fileQuadron = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_Break.dat";
	FILE* fileBreak = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_Combine.dat";
	FILE* fileCombine = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_neighborCell.dat";
	FILE* fileNeighborCell = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_coordinationNum.dat";
	FILE* fileCoordinationNum = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_coordinationNum0.dat";
	FILE* fileCoordinationNum0 = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_rattlerParticles.dat";
	FILE* fileRattler = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_cellMemberParticles.dat";
	FILE* fileMember = fopen(fileName.str().c_str(), "w");


	QuadronSystem::setLeesEdwardsBoundaryOffset(0.0);
	ifstream fileBoundary;
	fileBoundary.open(path + "output_boundary.dat");
	double xmin, xmax, ymin, ymax;
	fileBoundary >> xmin >> xmax >> ymin >> ymax;
	QuadronSystem::setPeriodicBoundarySize(xmin, xmax, ymin, ymax);
	QuadronSystem::QuadronAnalysis quadronAnalysis = QuadronSystem::QuadronAnalysis(path);
	quadronAnalysis.runAnalysis();

	quadronAnalysis.printNeighborCell2File(fileNeighborCell);
	quadronAnalysis.printCells2File(fileCell);
	quadronAnalysis.printParticle2File(fileParticle);
	quadronAnalysis.printQuadrons2File(fileQuadron);
	quadronAnalysis.printCoordinationNum2File(fileCoordinationNum);
	quadronAnalysis.printOriginalCoordinationNum2File(fileCoordinationNum0);
	quadronAnalysis.printRattlerParticle2File(fileRattler);
	quadronAnalysis.printCellMemberParticle2File(fileMember);


	fclose(fileCell);
	fclose(fileParticle);
	fclose(fileQuadron);
	fclose(fileBreak);
	fclose(fileCombine);
	fclose(fileNeighborCell);
	fclose(fileCoordinationNum);
	fclose(fileCoordinationNum0);
	fclose(fileRattler);
	fclose(fileMember);

	printf("finish\n");
	//std::cin.get();
}


int main()
{
	//std::vector<std::string> Alpha = {
	//	"alpha=1.0/",
	//	"alpha=1.1/",
	//	"alpha=1.2/",
	//	"alpha=1.3/",
	//	"alpha=1.4/",
	//	"alpha=1.5/",
	//	"alpha=1.6/",
	//	"alpha=1.7/",
	//	"alpha=1.8/",
	//	"alpha=1.9/",
	//	"alpha=2.0/",
	//	"alpha=2.5/",
	//	"alpha=3.0/",
	//	"alpha=4.0/" };
	std::vector<std::string> Alpha = {
		"alpha=1.01/",
		"alpha=1.05/"};
	std::vector<std::string> Mu = {
		"0/","1/","2/","3/","4/","5/"
	};

	for (const auto& alpha : Alpha) {
		for (const auto& mu : Mu) {
			std::stringstream path;
			path << alpha << mu << "Compression/Final/";
			oneCase(path.str().c_str());
		}
	}

	return 0;
}