#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <ctime>
#include <random>

#include <ios>     // std::left, std::right
#include <iomanip> // std::setw(int), std::setfill(char)

#include "src/CircularParticle.h"
#include "src/EllipticParticle.h"
#include "src/LeesEdwardsBoundaryWorld.h"
#include "src/QuadronAnalysis.h"

#pragma warning(disable : 4996)

int main()
{
	std::cout << "Simulation start..." << std::endl;
	std::cout << "Modification from my Mac..." << std::endl;

	FILE* filePara = fopen("parameters.dat", "r");
	if (!filePara) {
		printf("Can not open input parameter file...\n");
		return -1;
	}
	int numParticle;
	double aspectRatio;
	double friction;
	double pressure;
	double initialSolidFrac;
	int numThread, sizeChunk;
	fscanf(filePara, "%d", &numParticle);
	fscanf(filePara, "%lf", &aspectRatio);
	fscanf(filePara, "%lf", &friction);
	fscanf(filePara, "%lf", &pressure);
	fscanf(filePara, "%lf", &initialSolidFrac);
	fscanf(filePara, "%d", &numThread);
	fscanf(filePara, "%d", &sizeChunk);
	printf("number of particle: %d\n", numParticle);
	printf("friction: %lf\n", friction);
	printf("pressure: %lf\n", pressure);
	printf("OMP parameters:\nnumber of thread: %d\nsize of chunk: %d\n", numThread, sizeChunk);


	double dt = 1e-4;
	const double density = 10.0;
	const double kn = 1e7;
	const double kt = 1e7;
	const double resCoeff = 0.01;
	const double friCoeff = friction;
	const double xmin = -100.0;
	const double xmax = 100.0;
	const double ymin = -100.0;
	const double ymax = 100.0;

	const double radiusScaleFactor = sqrt(aspectRatio);

	std::default_random_engine engine(clock());
	std::uniform_real_distribution<double> real_r(2.0 / 3.0, 4.0 / 3.0);
	std::uniform_real_distribution<double> real_x(xmin, xmax);
	std::uniform_real_distribution<double> real_y(ymin, ymax);
	std::uniform_real_distribution<double> real_t(0.0, M_PI * 2.0);

	// RANDOM GENERATION
	// Generate particles, contact model, and world.
	ContactModel contactModel(kn, kt, 0.95, 0.0);
	vector<CircularParticle> circularParticles(numParticle);
	vector<BaseParticle*> particleHandlers(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		circularParticles[i] = CircularParticle(i, real_r(engine) * radiusScaleFactor,
			vector2d(real_x(engine), real_y(engine)), vector2d(0.0, 0.0),
			real_t(engine), 0.0, &contactModel);
		circularParticles[i].calculateParticleProperties(density);
		particleHandlers[i] = &circularParticles[i];
	}
	circularParticles[0].setPos(vector2d(0.0, 0.0));
	circularParticles[0].appendSetting(SKIP_TIME_INTEGRAL);
	PeriodicBoundaryWorld world_rg(dt, particleHandlers, &contactModel);
	world_rg.set_OMP_NumThreadsAndSizeChunk(numThread, sizeChunk);
	world_rg.setWorldBoundary(xmin, xmax, ymin, ymax);
	// Re-calculate boundary size to satisfy given solid fraction.
	world_rg.scaleBoundarySizeToGivenSolidFraction(initialSolidFrac);

	// Simulation start...
	string path = "RandomGeneration\\";
	world_rg.prepare(path);
	for (int iStep = 0; iStep < 5e4; ++iStep) {
		world_rg.takeTimeStep();
		world_rg.resetVelocityOfRattlers();
		if (fmod(iStep, 1000) == 0) {
			world_rg.print2Screen_worldState(iStep);
			world_rg.writeParticleTimeHistory2Files();
			world_rg.flushAllFiles();
		}
		if (world_rg.countActucalContacts() < 5) break;
	}
	world_rg.closeParticleTimeHistoryFiles();
	printf("Finish...\n");




	// Generate elliptic partiles
	vector<EllipticParticle> ellipticParticles(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		ellipticParticles[i] = EllipticParticle(circularParticles[i], aspectRatio, &contactModel);
		ellipticParticles[i].calculateParticleProperties(density);
		particleHandlers[i] = &ellipticParticles[i];
	}
	ellipticParticles[0].appendSetting(SKIP_TIME_INTEGRAL);

	// Control parameters...
	double ConstantStress = pressure;
	// Generate particles and world
	PeriodicControlWorld world_le(dt, particleHandlers, &contactModel);
	world_le.set_OMP_NumThreadsAndSizeChunk(numThread, sizeChunk);
	world_le.setWorldBoundary(world_rg.getXmin(), world_rg.getXmax(), world_rg.getYmin(), world_rg.getYmax());
	world_le.updateTotalParticleVolume();
	world_le.setVirtualBoundaryMass();
	// Simulation start...
	// Compression
	path = "Compression\\";
	world_le.prepare(path);
	world_le.closeParticleTimeHistoryFiles();
	contactModel.setFriction(friCoeff);
	contactModel.setRestitutionCoeff(0.01);
	contactModel.recalculateFactor();
	int iStep = 0;
	while (1) {


		world_le.scaleVelocityOfRattlers(0.9999);
		world_le.modifyParticlePosition();
		world_le.findPossibleContacts();
		world_le.updateContacts();
		world_le.collectForceAndTorque();
		world_le.takeTimeIntegral();
		world_le.updateTotalStress();
		//world_le.updatePeriodicBoundary_stressControlUniform(ConstantStress);
		world_le.updatePeriodicBoundary_stressControlX(ConstantStress);
		world_le.updatePeriodicBoundary_stressControlY(ConstantStress);


		if (fmod(iStep, 10000) == 0) {
			char path0[50];
			sprintf(path0, "%s%010d", path.c_str(), iStep);
			world_le.prepare2(path0);
			world_le.print2Screen_worldState(iStep);
			world_le.writeParticleTimeHistory2Files();
			world_le.flushAllFiles();
			world_le.closeParticleTimeHistoryFiles();
		}

		if ((world_le.getKineticEnergyPerParticle() < 1e-12 && iStep > 1e6)) break;

		iStep++;
	}






	// Final step output
	path = "Compression\\Final\\";
	world_le.prepare(path);
	world_le.writeParticleTimeHistory2Files();
	world_le.flushAllFiles();
	world_le.closeParticleTimeHistoryFiles();


	printf("==================================================\n");
	printf("================ SIMULATION FINISH ===============\n");
	printf("==================================================\n");

	return 0;
}
