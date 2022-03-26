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
	const double resCoeff = 1.0;
	const double friCoeff = friction;
	const double xmin = -100.0;
	const double xmax = 100.0;
	const double ymin = -100.0;
	const double ymax = 100.0;

	const double radiusScaleFactor = sqrt(aspectRatio);

	std::default_random_engine engine(0);
	std::uniform_real_distribution<double> real_r(2.0 / 3.0, 4.0 / 3.0);
	std::uniform_real_distribution<double> real_x(xmin, xmax);
	std::uniform_real_distribution<double> real_y(ymin, ymax);
	std::uniform_real_distribution<double> real_t(0.0, M_PI * 2.0);
	std::uniform_real_distribution<double> real_a(1.01, aspectRatio);

	// RANDOM GENERATION
	// Generate particles, contact model, and world.
	ContactModel contactModel(kn, kt, 0.95, 0.0);
	vector<CircularParticle> circularParticles(numParticle);
	vector<BaseParticle*> particleHandlers(numParticle);
	vector<double> alpha(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		alpha[i] = real_a(engine);
		circularParticles[i] = CircularParticle(i, real_r(engine) * radiusScaleFactor,
			vector2d(real_x(engine), real_y(engine)), vector2d(0.0, 0.0),
			real_t(engine), 0.0, &contactModel);
		circularParticles[i].calculateParticleProperties(density);
		particleHandlers[i] = &circularParticles[i];
	}
	//circularParticles[0].setPos(vector2d(0.0, 0.0));
	//circularParticles[0].appendSetting(SKIP_TIME_INTEGRAL);
	PeriodicBoundaryWorld simRandGen(dt, particleHandlers, &contactModel);
	simRandGen.set_OMP_NumThreadsAndSizeChunk(numThread, sizeChunk);
	simRandGen.setWorldBoundary(xmin, xmax, ymin, ymax);
	// Re-calculate boundary size to satisfy given solid fraction.
	simRandGen.scaleBoundarySizeToGivenSolidFraction(initialSolidFrac);

	// Simulation start...
	string path = "RandomGeneration\\";
	simRandGen.prepare(path);
	for (int iStep = 0; iStep < 5e4; ++iStep) {
		simRandGen.takeTimeStep();
		simRandGen.resetVelocityOfRattlers();
		if (fmod(iStep, 1000) == 0) {
			simRandGen.print2Screen_worldState(iStep);
			simRandGen.writeParticleTimeHistory2Files();
			simRandGen.flushAllFiles();
		}
		if (simRandGen.countActucalContacts() < 5) break;
	}
	simRandGen.closeParticleTimeHistoryFiles();
	printf("Finish...\n");




	// Generate elliptic partiles
	vector<EllipticParticle> ellipticParticles(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		ellipticParticles[i] = EllipticParticle(circularParticles[i], aspectRatio, &contactModel);
		ellipticParticles[i].calculateParticleProperties(density);
		particleHandlers[i] = &ellipticParticles[i];
	}
	//ellipticParticles[0].appendSetting(SKIP_TIME_INTEGRAL);

	// Control parameters...
	double ConstantStress = pressure;
	// Generate particles and world
	PeriodicControlWorld simIsoComp(dt, particleHandlers, &contactModel);
	simIsoComp.set_OMP_NumThreadsAndSizeChunk(numThread, sizeChunk);
	simIsoComp.setWorldBoundary(simRandGen.getXmin(), simRandGen.getXmax(), simRandGen.getYmin(), simRandGen.getYmax());
	simIsoComp.updateTotalParticleVolume();
	simIsoComp.setVirtualBoundaryMass();
	// Simulation start...
	// Compression
	path = "Compression\\";
	simIsoComp.prepare(path);
	contactModel.setFriction(friCoeff);
	contactModel.setRestitutionCoeff(0.01);
	contactModel.recalculateFactor();

	FILE* fileTime = fopen("Compression\\output_computationalTime.dat", "w");
	clock_t timer = clock();
	int iStep = 0;
	int count = 0;
	while (1) {

		simIsoComp.scaleVelocityOfRattlers(0.9);
		simIsoComp.modifyParticlePosition();
		simIsoComp.findPossibleContacts();
		simIsoComp.updateContacts();
		simIsoComp.collectForceAndTorque();
		simIsoComp.takeTimeIntegral();
		simIsoComp.updateTotalStress();
		if (simIsoComp.countActucalContacts() == 0 || simIsoComp.getKineticEnergyPerNonRattlerParticle() < 1e-8) {
			simIsoComp.updatePeriodicBoundary_strainControl(-1e-2, -1e-2);
			count++;
		}


		if (fmod(iStep, 10000) == 0) {
			simIsoComp.print2Screen_worldState(iStep);
		}

		if (count == 10) {
			simIsoComp.writeParticleTimeHistory2Files();
			simIsoComp.flushAllFiles();

			timer = clock() - timer;
			fprintf(fileTime, "%lf\n", static_cast<double>(timer) / CLOCKS_PER_SEC);
			timer = clock();
			fflush(fileTime);

			count = 0;
		}

		iStep++;
	}
	
	simIsoComp.closeParticleTimeHistoryFiles();

	path = "Compression\\Final\\";
	simIsoComp.prepare(path);
	simIsoComp.writeParticleTimeHistory2Files();
	simIsoComp.flushAllFiles();
	simIsoComp.closeParticleTimeHistoryFiles();

	fclose(fileTime);


	printf("==================================================\n");
	printf("================ SIMULATION FINISH ===============\n");
	printf("==================================================\n");

	return 0;
}
