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


	printf("finish\n");
	std::cin.get();
}


int main()
{

	oneCase("");

	return 0;
}