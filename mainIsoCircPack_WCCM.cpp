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



void Simulation(const std::string path0, 
	const int numParticle,
	const double aspectRatio,
	const double friction,
	const double pressure,
	const int numThread,
	const int sizeChunk)
{
	std::cout << "Simulation start..." << std::endl;
	//std::cout << "Modification from my Mac..." << std::endl;

	

	double dt = 1e-4;
	double initialSolidFrac = 0.8;
	const double density = 10.0;
	const double kn = 1e7;
	const double kt = 1e7;
	const double resCoeff = 0.01;
	const double friCoeff = friction;
	const double xmin = -100.0;
	const double xmax = 100.0;
	const double ymin = -100.0;
	const double ymax = 100.0;

	std::default_random_engine engine(0);
	std::uniform_real_distribution<double> real_r(2.0 / 3.0, 4.0 / 3.0);
	std::uniform_real_distribution<double> real_x(xmin, xmax);
	std::uniform_real_distribution<double> real_y(ymin, ymax);
	std::uniform_real_distribution<double> real_t(0.0, M_PI * 2.0);

	// RANDOM GENERATION
	// Generate particles, contact model, and world.
	ContactModel contactModel(kn, kt, 0.95, 0.0);
	vector<CircularParticle> circularParticles(numParticle);
	vector<BaseParticle*> particleHandlers(numParticle);
	vector<double> alpha(numParticle);
	for (int i = 0; i < numParticle; ++i) {
		circularParticles[i] = CircularParticle(i, real_r(engine),
			vector2d(real_x(engine), real_y(engine)), vector2d(0.0, 0.0),
			real_t(engine), 0.0, &contactModel);
		circularParticles[i].calculateParticleProperties(density);
		particleHandlers[i] = &circularParticles[i];
	}
	//circularParticles[0].setPos(vector2d(0.0, 0.0));
	//circularParticles[0].appendSetting(SKIP_TIME_INTEGRAL);
	PeriodicBoundaryWorld world_rg(dt, particleHandlers, &contactModel);
	world_rg.set_OMP_NumThreadsAndSizeChunk(numThread, sizeChunk);
	world_rg.setWorldBoundary(xmin, xmax, ymin, ymax);
	// Re-calculate boundary size to satisfy given solid fraction.
	world_rg.scaleBoundarySizeToGivenSolidFraction(initialSolidFrac);

	// Simulation start...
	string path = path0 + "RandomGeneration\\";
	world_rg.prepare(path);
	for (int iStep = 0; iStep < 5e4; ++iStep) {
		world_rg.takeTimeStep();
		world_rg.resetVelocityOfRattlers();
		if (fmod(iStep, 100) == 0) {
			world_rg.print2Screen_worldState(iStep);
			world_rg.writeParticleTimeHistory2Files();
			world_rg.flushAllFiles();
		}
		if (world_rg.countActucalContacts() < 20) break;
	}
	world_rg.closeParticleTimeHistoryFiles();
	printf("Finish...\n");

	FILE* fileMu = fopen((path0 + "output_mu.dat").c_str(), "w");
	fprintf(fileMu, "%lf", friction);
	fclose(fileMu);



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
	path = path0 + "Compression\\";
	world_le.prepare(path);
	contactModel.setFriction(friCoeff);
	contactModel.setRestitutionCoeff(resCoeff);
	BaseParticle::_GlobalDamping_ = 0.1;
	contactModel.recalculateFactor();

	FILE* fileTime = fopen((path + "output_computationalTime.dat").c_str(), "w");
	clock_t timer = clock();
	int iStep = 0;
	while (1) {

		world_le.scaleVelocityOfRattlers(0.99);
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
			world_le.print2Screen_worldState(iStep);
			world_le.writeParticleTimeHistory2Files();
			world_le.flushAllFiles();

			timer = clock() - timer;
			fprintf(fileTime, "%lf\n", static_cast<double>(timer) / CLOCKS_PER_SEC);
			timer = clock();
			fflush(fileTime);
		}

		if ((world_le.getKineticEnergyPerParticle() < 1e-10 && iStep > 1e5)) break;

		iStep++;
	}
	
	world_le.closeParticleTimeHistoryFiles();

	path = path0 + "Compression\\Final\\";
	world_le.prepare(path);
	world_le.writeParticleTimeHistory2Files();
	world_le.flushAllFiles();
	world_le.closeParticleTimeHistoryFiles();

	fclose(fileTime);


	printf("==================================================\n");
	printf("================ SIMULATION FINISH ===============\n");
	printf("==================================================\n");

}




int main()
{
	FILE* filePara = fopen("parameters.dat", "r");
	if (!filePara) {
		printf("Can not open input parameter file...\n");
		return -1;
	}
	int numParticle;
	double aspectRatio;
	double pressure;
	int numThread, sizeChunk;
	fscanf(filePara, "%d", &numParticle);
	fscanf(filePara, "%lf", &aspectRatio);
	fscanf(filePara, "%lf", &pressure);
	fscanf(filePara, "%d", &numThread);
	fscanf(filePara, "%d", &sizeChunk);
	printf("number of particle: %d\n", numParticle);
	printf("pressure: %lf\n", pressure);
	printf("OMP parameters:\nnumber of thread: %d\nsize of chunk: %d\n", numThread, sizeChunk);


	int caseNum;
	std::vector<double> frictions;
	double friction;
	fscanf(filePara, "%d", &caseNum);
	for (int i = 0; i < caseNum; ++i) {
		fscanf(filePara, "%lf", &friction);
		frictions.push_back(friction);
	}

	std::cout << "Read para finish..." << std::endl;
	std::cout << frictions.size() << std::endl;

	for (int i = 0; i < frictions.size(); ++i) {
		std::cout << i << std::endl;
		std::stringstream path;
		path << i << "\\";
		std::cout << path.str() << std::endl;
		Simulation(path.str(), numParticle, aspectRatio, frictions[i], pressure, numThread, sizeChunk);
	}


	return 0;
}