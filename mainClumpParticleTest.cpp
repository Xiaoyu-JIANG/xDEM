#include <ctime>
#include <random>
#include <sstream>
#include <fstream>

#include "src/CircularParticle.h"
#include "src/EllipticParticle.h"
#include "src/GluedParticle.h"
#include "src/PeriodicControlWorld.h"

#pragma warning(disable : 4996)

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::stringstream;

int main() {

	// create output folders
	string path = "debugOutput\\";
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());

	// Define simulation parameters
	int numParticle = 1000;
	const double dt = 1e-4;
	const double density = 1.0;
	const double kn = 1e7;
	const double kt = 1e7;
	const double resCoeff = 0.8;
	const double friCoeff = 0.01;
	const double xmin = -20.0;
	const double xmax = 20.0;
	const double ymin = -20.0;
	const double ymax = 20.0;
	const double initialSolidFrac = 0.80;

	// Random number generations used later
	std::default_random_engine engine(clock());
	std::uniform_int_distribution<int> int_num(1, 4);
	std::uniform_real_distribution<double> real_r(1.0, 2.0);
	std::uniform_real_distribution<double> real_k(1.1, 4.0);
	std::uniform_real_distribution<double> real_theta(0.0, 2.0 * M_PI);
	std::uniform_real_distribution<double> real_l(-0.25, 0.25);
	std::uniform_real_distribution<double> real_x(xmin, xmax);
	std::uniform_real_distribution<double> real_y(ymin, ymax);

	// Create contact model objects
	ContactModel contactModel(kn, kt, resCoeff, friCoeff);
	// Create circular particles in initial random generation
	vector<CircularParticle> circularParticles0(numParticle);
	vector<BaseParticle*> particleHandlers0(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		circularParticles0[i] = CircularParticle(i, real_r(engine),
			vector2d(real_x(engine), real_y(engine)), vector2d(0.0, 0.0),
			real_theta(engine), 0.0, &contactModel);
		circularParticles0[i].calculateParticleProperties(density);
		particleHandlers0[i] = &circularParticles0[i];
	}
	// Link particles and contact model to world
	PeriodicBoundaryWorld world_randomGeneration(dt, particleHandlers0, &contactModel);
	world_randomGeneration.setWorldBoundary(xmin, xmax, ymin, ymax);
	world_randomGeneration.set_OMP_NumThreadsAndSizeChunk(6, 100);
	// Re-calculate boundary size to satisfy given solid fraction.
	world_randomGeneration.scaleBoundarySizeToGivenSolidFraction(initialSolidFrac);



	// Simulation start...
	path = "debugOutput\\RandomGeneration\\";
	world_randomGeneration.prepare(path);
	for (int iStep = 0; iStep < 3e4; ++iStep) {

		world_randomGeneration.modifyParticlePosition();
		world_randomGeneration.findPossibleContacts();
		world_randomGeneration.updateContacts();
		world_randomGeneration.collectForceAndTorque();
		world_randomGeneration.takeTimeIntegral();

		if (fmod(iStep, 1000) == 0) {
			world_randomGeneration.writeParticleTimeHistory2Files();
			world_randomGeneration.flushAllFiles();
			world_randomGeneration.print2Screen_worldState(iStep);
		}
	}
	world_randomGeneration.closeParticleTimeHistoryFiles();





	// create elliptic particles for glued particles
	vector<vector<EllipticParticle>> ellipticParticles(numParticle);
	vector<GluedParticle> gluedParticles(numParticle);
	vector<BaseParticle*> particleHandlers;
	int countID = 0;
	for (int i = 0; i < numParticle; ++i) {
		int numSubParticle = int_num(engine);
		ellipticParticles[i].resize(numSubParticle);
		if (numSubParticle == 1) {
			ellipticParticles[i][0] = EllipticParticle(countID,
				circularParticles0[i].getOuterRadius(), real_k(engine),
				circularParticles0[i].getPos(),
				vector2d(0.0, 0.0),
				real_theta(engine),
				0.0,
				&contactModel);
			particleHandlers.push_back(&ellipticParticles[i][0]);
			countID++;
			ellipticParticles[i][0].calculateParticleProperties(density);
		}
		else {
			// Generate main particle
			gluedParticles[i] = GluedParticle(countID, vector2d(0.0, 0.0), 0.0);
			particleHandlers.push_back(&gluedParticles[i]);
			countID++;
			// glued particle itself is not involed in contact detection
			gluedParticles[i].appendSetting(SETTING::SKIP_CONTACT_DETECTION);
			// Generate sub particles
			double tempTheta = real_theta(engine);
			for (int j = 0; j < numSubParticle; ++j) {
				ellipticParticles[i][j] = EllipticParticle(countID,
					circularParticles0[i].getOuterRadius(), 3.0,
					circularParticles0[i].getPos() + vector2d(real_l(engine), real_l(engine)),
					vector2d(0.0, 0.0),
					tempTheta + M_PI * j / numSubParticle,
					0.0,
					&contactModel);
				particleHandlers.push_back(&ellipticParticles[i][j]);
				countID++;
				ellipticParticles[i][j].setUpLevel(gluedParticles[i].getID());
				// sub-particles are not involed in time integral, their states (position, velocity...) are modified through their main-particle
				ellipticParticles[i][j].appendSetting(SETTING::SKIP_TIME_INTEGRAL);
				ellipticParticles[i][j].appendSetting(SETTING::SKIP_MACRO_SUM);
				gluedParticles[i].pushBackSubParticle(&ellipticParticles[i][j]);
			}
			// Calculate glued particle properties
			gluedParticles[i].calculateParticleProperties(density);
		}

	}
	// Link new particle handlers to new world
	PeriodicControlWorld world_uniformCompression(dt, particleHandlers, &contactModel);
	world_uniformCompression.set_OMP_NumThreadsAndSizeChunk(8, 100);
	world_uniformCompression.setWorldBoundary(world_randomGeneration.getXmin(), world_randomGeneration.getXmax(), world_randomGeneration.getYmin(), world_randomGeneration.getYmax());
	world_uniformCompression.updateTotalParticleVolume();
	world_uniformCompression.setVirtualBoundaryMass();
	// Control parameters...
	double currentSolidFrac = world_uniformCompression.getSolidFraction();
	// Simulation start...
	path = "debugOutput\\Compression\\";
	world_uniformCompression.prepare(path);


	for (int iStep = 0; iStep < 1e5; ++iStep) {


		world_uniformCompression.modifyParticlePosition();
		world_uniformCompression.findPossibleContacts();
		world_uniformCompression.updateContacts();
		world_uniformCompression.collectForceAndTorque();
		world_uniformCompression.takeTimeIntegral();
		world_uniformCompression.resetVelocityOfRattlers();
		world_uniformCompression.updateTotalStress();
		world_uniformCompression.updatePeriodicBoundary_stressControlUniform(1e3);

		if (fmod(iStep, 1000) == 0) {
			world_uniformCompression.writeParticleTimeHistory2Files();
			world_uniformCompression.flushAllFiles();
			world_uniformCompression.print2Screen_worldState(iStep);
		}
	}
	world_uniformCompression.closeParticleTimeHistoryFiles();

	return 0;
}
