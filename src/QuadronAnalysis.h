#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>

#include "BaseParticle.h"
#include "MathMatrix2d.h"

using std::string;
using std::to_string;
using std::pair;
using std::vector;
using std::unordered_map;
using std::ifstream;


namespace QuadronSystem {

	class Polygon;
	class Particle;
	class Contact;
	class Quadron;
	class Cell;
	class QuadronAnalysis;

	// For periodic boundary modifications.
	extern double _xmin, _xmax, _ymin, _ymax, _xlen, _ylen, _xlenHalf, _ylenHalf;
	// For Lees-Edwards boundary modifications.
	extern double _offset;
	// Set parameters of periodic boundary modification used in QuadronSystem scope...
	void setPeriodicBoundarySize(const double& xmin, const double& xmax, const double& ymin, const double& ymax);
	// Set parameters of Lees-Edwards boundary modification used in QuadronSystem scope...
	void setLeesEdwardsBoundaryOffset(const double& offset);
	// Modify temporal position due to boundary conditions...
	void modifyPosition(double& tempX, double& tempY, const double& x, const double& y);
	void print();
	// Do quadron analysis for single step...
	void quadronAnalysisInThisStep(const string path, vector<BaseParticle*>& particleHandlers);

	extern int NUM_THREADS;
	extern int SIZE_CHUNK;

	extern unordered_map<string, pair<int, int>> oldContactCellList;
	extern unordered_map<string, pair<int, int>> newContactCellList;
	extern unordered_map<string, Polygon> oldCellList;
	extern unordered_map<string, Polygon> newCellList;
	extern unordered_map<string, pair<int, int>> cellSurviveTime;

	// ==========================================
	// [0] Polygon object
	// ==========================================
	class Polygon
	{
	public:
		vector<double> _x;
		vector<double> _y;
		double _cx = 0.0;
		double _cy = 0.0;
		double _area = 0.0;

	public:
		Polygon();
		void pushBackVertex(const double& x, const double& y);
		void calculateProperties();
		int pointInsidePolygon(const double& x, const double& y);
		void print2File(FILE* file);
	};

	// ==========================================
	// [1] Particle object
	// ==========================================
	class Particle
	{
	public:
		int _id = -1;
		bool _isActive = true;
		// Position of real particle.
		double _x = 0.0;
		double _y = 0.0;
		double _volume = 0.0;
		// Neighbouring particle and contact list.
		vector<QuadronSystem::Particle*> _pParticles;
		vector<QuadronSystem::Contact*> _pContacts;
		// Sorted neighbouring particle and contact list.
		// They are not complete, inactivated particles and contacts are not included.
		// "_link" is used to specify if the path to next particle is accessibale.
		vector<QuadronSystem::Particle*> _pParticlesSorted;
		vector<QuadronSystem::Contact*> _pContactsSorted;
		vector<int> _link;
		// Polygon shape particle
		QuadronSystem::Polygon _polygon;
		vector<QuadronSystem::Quadron*> _pQuadrons;
		matrix2d _stress;
		double _totalArea = 0.0;

	public:
		Particle();
		Particle(const int& id, const double& x, const double& y, const double& vol);
		void pushBackSurroundingParticle(QuadronSystem::Particle* P);
		void pushBackSurroundingContact(QuadronSystem::Contact* c);
		void setAsInactiveParticle();
		void setAllSurroundingContactsAsInactive();
		void rearrangeContactsClockwisely();
		void constructPolygon();
		void calculateStress();
	};


	// ==========================================
	// [2] Contact object
	// ==========================================
	class Contact
	{
	public:
		int _id = -1;
		bool _isActive = true;
		// Position of contact point.
		double _x = 0.0;
		double _y = 0.0;
		Particle* _pParticle1 = nullptr;
		Particle* _pParticle2 = nullptr;
		vector<Cell*> _pCell;
		vector2d _force;

	public:
		Contact(const int& id, const double& x, const double& y,
			QuadronSystem::Particle* p1, QuadronSystem::Particle* p2,
			const vector2d& force);
		void setAsInactiveContact();
	};


	// ==========================================
	// [3] Quadron object
	// ==========================================
	class Quadron
	{
	public:
		int _id = -1;
		int _cellOrder = 0;
		Particle* _pParticle = nullptr;
		Cell* _pCell = nullptr;
		Contact* _pContact1 = nullptr;
		Contact* _pContact2 = nullptr;
		Polygon _polygon;

	public:
		Quadron(const int id,
			QuadronSystem::Particle* p, QuadronSystem::Cell* c,
			QuadronSystem::Contact* c1, QuadronSystem::Contact* c2);
		void constructPolygon();
		void setPointerToParticlesAndCells();
	};


	// ==========================================
	// [4] Cell object
	// ==========================================
	class Cell
	{
	public:
		int _id = -1;
		int _cellOrder = 0;
		// List of particles and quadrons constructing this cell.
		vector<QuadronSystem::Particle*> _pParticles;
		vector<QuadronSystem::Contact*> _pContacts;
		vector<QuadronSystem::Quadron*> _pQuadrons;
		QuadronSystem::Polygon _polygon;
		matrix2d _stress;
		double _totalArea = 0.0;

	public:
		Cell(QuadronSystem::Particle* particle, QuadronSystem::Contact* contact);
		void pushBackNewParticle(QuadronSystem::Particle* particle);
		void pushBackNewContact(QuadronSystem::Contact* contact);
		void pushBackNewEdge(QuadronSystem::Particle* particle, QuadronSystem::Contact* contact);
		void setPointerToContacts();
		void constructPolygon();
		void calculateStress();
	};


	// ==========================================
	// [5] Analysis class
	// ==========================================
	class QuadronAnalysis
	{
	private:
		vector<QuadronSystem::Particle> particles;
		vector<QuadronSystem::Contact> contacts;
		vector<QuadronSystem::Quadron> quadrons;
		vector<QuadronSystem::Cell> cells;

	public:
		// [0] Constructors... input: DEM particles
		QuadronAnalysis(const string path);
		QuadronAnalysis(const vector<BaseParticle*>& p);
		QuadronAnalysis(vector<double>& vol, ifstream& filePos, ifstream& fileContact);

		void runAnalysis();

		void makeCellList();
		void compareCellLists(FILE* fileB, FILE* fileC);
		void checkCellSurviveTime(FILE* file);

		/*DO NOT USE THIS FUNCTION, USE makeCellList...*/
		void makeContactCellList();
		/*DO NOT USE THIS FUNCTION, USE compareCellList*/
		void compareContactCellLists(FILE* fileB, FILE* fileC, const bool& isPrint2File);

		double calculateMeanCoordinationNumber();
		double calculateMeanOriginalCoordinationNumber();
		double calculateRattlerFreeSolidFraction();
		double calculateSolidFraction();
		vector<double> calculateDetailCoordinationNumber();
		vector<double> calculateDetailOriginalCoordinationNumber();

		// print information to file...
		void printQuadrons2File(const std::string& path);
		void printQuadrons2File(FILE* file);
		void printCells2File(const std::string& path);
		void printCells2File(FILE* file);
		void printParticle2File(FILE* file);
		void printNeighborCell2File(FILE* file);
		void printCoordinationNum2File(FILE* file);
		void printOriginalCoordinationNum2File(FILE* file);
		void printSolidFraction2File(FILE* file);
		void printRattlerParticle2File(FILE* file);
		void printCellMemberParticle2File(FILE* file);

		// Check quadron area...
		double checkQudronArea();
		// Check particle area...
		double checkParticleArea();
		// Check cell area...
		double checkCellArea();

		// get method
		int getCellNum();

	private:
		// [1] Pre-inactivate rattler particles with no contact...
		void inactivateRattlerParticles();
		// [2] Remove the open chain structures...
		// Removement of particle and contact will influence its neighbour structure
		// Thus each time only one particle and one contact are removed
		void removeOpenChainParticles();
		// [3] Reordering the neighbouring particle and contact list...
		// From the point of view of particle center
		// its neighbouring particles and contacts are located in a clockwise order
		void rearrangeContactsClockwisely();
		// ===============================================
		// [4] MAIN ALGORITHM OF CELL(QUADRON) CONSTRUCTION...
		// ===============================================
		void constructCellsFromParticlesAndContacts();
		void constructCellsFromParticlesAndContacts_anotherVersion();
		void giveCellID();
		// [5] Construct quadrons from cells...
		void constructQuadronsFromCells();
		// [6] Calculate particle and cell stress...
		void calculateParticleStress();
		void calculateCellStress();
		
		// Check all particles if it have only one contact...
		// Return the pointer of the firstly found particle
		// This function is used in "removeOpenChainParticles"
		QuadronSystem::Particle* findParticleWithOnlyOneContact();
		QuadronSystem::Particle* findParticleWithLessThanThreeContacts();
		
		
	};

}

