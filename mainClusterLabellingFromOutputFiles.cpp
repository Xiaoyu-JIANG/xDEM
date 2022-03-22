#include "src/ClusterLabelling.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


using namespace std;

int main()
{
	cout << "Hello, cluster labelling..." << endl;
	const string path = "";
	ClusterSystem::ClusterAnalysis pack = ClusterSystem::ClusterAnalysis(path, 1, M_PI / 180);
	pack.runAnalysis();
	pack.printClusters2File(path);


	return 0;
}