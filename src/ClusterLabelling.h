#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <limits>

#include "BaseParticle.h"
#include "MathMatrix2d.h"

using std::string;
using std::to_string;
using std::pair;
using std::vector;
using std::unordered_map;
using std::ifstream;
using std::stringstream;

namespace ClusterSystem
{
	class Particle;
	class Contact;
	class Cluster;
	class ClusterAnalysis;


	class Particle
	{
	public:
		int _id = -1;
		double _theta;
		
	public:
		Particle(const int id, const double theta);
	};

	class Contact
	{
	public:
		double _dTheta;
		double _normalForce;
		double _shearForce;
		double _totalForce;
		Particle* _particle1 = nullptr;
		Particle* _particle2 = nullptr;

	public:
		Contact(Particle* p1, Particle* p2, double fn, double fs, double f);
		void calculateRelativeAngle();
	};

	class Cluster
	{
	public:
		int _id;
		Cluster* _pnode = nullptr;

		Cluster() :_id(-1) {}
		Cluster(int id);

		int findHeadParentID();
		void linekTo(Cluster& parent);
	};

	class ClusterAnalysis
	{
	public:
		vector<Particle> _particles;
		vector<Contact> _contacts;
		vector<double> _dTheta0;
		vector<vector<int>> _nclus;
		vector<vector<int>> _mclus;
		vector<double> _correlation;

	public:
		ClusterAnalysis(const vector<BaseParticle*>& p, const double accuracy);
		ClusterAnalysis(const string path, const int iStep, const double accuracy);

		void runAnalysis();
		void printClusters2File(string path);

	private:
		void calculateRelativeAngle();
	};
}