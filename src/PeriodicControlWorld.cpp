#include <omp.h>
#include <sstream>
#include <ctime>
#include "PeriodicControlWorld.h"

using std::stringstream;

//void PeriodicControlWorld::setControlMethod(
//	const ControlMethod& controlMethodX, const double& valueX, 
//	const ControlMethod& controlMethodY, const double& valueY)
//{
//	_controlMethodX = controlMethodX;
//	_controlMethodY = controlMethodY;
//	printf("\n[MESSAGE]\n");
//	printf("Periodic boundary control method changed...\n");
//	switch (_controlMethodX) {
//	case ControlMethod::Null:
//		printf("[ERROR]___X___Null, Please set X periodic boundary control method properly...\n");
//		exit(-1);
//	case ControlMethod::Fixed:
//		printf("___X___Fixed\n");
//		break;
//	case ControlMethod::Strain:
//		_strainRateX = valueX;
//		printf("___X___Strain control: %le\n", _strainRateX);
//		break;
//	case ControlMethod::Stress:
//		_stressX = valueX;
//		printf("___X___Stress control: %le\n", _stressX);
//		break;
//	default:
//		printf("[ERROR] ___X___No match control method...\n");
//		exit(-1);
//	}
//	switch (_controlMethodY) {
//	case ControlMethod::Null:
//		printf("[ERROR]___X___Null, Please set Y periodic boundary control method properly...\n");
//		exit(-1);
//	case ControlMethod::Fixed:
//		printf("___Y___Fixed\n");
//		break;
//	case ControlMethod::Strain:
//		_strainRateY = valueY;
//		printf("___Y___Strain control: %le\n", _strainRateY);
//		break;
//	case ControlMethod::Stress:
//		_stressY = valueY;
//		printf("___Y___Stress control: %le\n", _stressY);
//		break;
//	default:
//		printf("[ERROR] ___Y___No match control method...\n");
//		exit(-1);
//	}
//}

void PeriodicControlWorld::setVirtualBoundaryMass()
{
	_mass = 0.0;
	for (auto& iter : _particleHandlers) {
		if (iter->getSettings() & SETTING::SKIP_MACRO_SUM) continue;
		_mass += iter->getMass();
	}
	_mass /= sqrt(double(_particleHandlers.size()));
}

void PeriodicControlWorld::setVirtualBoundaryMass(const double& mass)
{
	_mass = mass;
}

void PeriodicControlWorld::setMaxStrainRate(const double maxStrainRate)
{
	_maxStrainRate = maxStrainRate;
}

void PeriodicControlWorld::print2Screen_worldState(const int& iStep)
{
	static clock_t timerTotal = clock();
	static clock_t timer = clock();
	timer = clock() - timer;
	printf("Current step = %d, Time used from lase output steps = %lf, Total time used = %lf\n", iStep,
		(static_cast<double>(timer) / CLOCKS_PER_SEC),
		((static_cast<double>(clock()) - timerTotal) / CLOCKS_PER_SEC));
	timer = clock();
	printf("\tCurrent solid fraction = %lf\n", getSolidFraction());
	printf("\tPossible contact number = %d\n", countPossibleContacts());
	printf("\tActual contact number = %d\n", countActucalContacts());
	printf("\tSystem kinetic energy per particle = %le\n", getKineticEnergyPerParticle());
	printf("\tSystem kinetic energy per non-rattler particle = %le\n", getKineticEnergyPerNonRattlerParticle());
	printf("\tSystem hydro stress = %le\n", getStress().trace() / 2.0);
}

void PeriodicControlWorld::prepare(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());
	checkTimeIncrement();
	writeParticleInformation2Files(path);
	createParticleTimeHistoryFiles(path);
}

void PeriodicControlWorld::prepare2(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());
	writeParticleInformation2Files(path);
	createParticleTimeHistoryFiles(path);
}

void PeriodicControlWorld::takeTimeStep()
{
	modifyParticlePosition();
	findPossibleContacts();
	updateContacts();
	collectForceAndTorque();
	takeTimeIntegral();
	resetVelocityOfRattlers();
}

void PeriodicControlWorld::updatePeriodicBoundary_strainControl(const double& strainRateX, const double& strainRateY)
{
	double scaleX = 1.0 + strainRateX * _xlen * _dt;
	double scaleY = 1.0 + strainRateY * _ylen * _dt;
	_xmin *= scaleX;
	_xmax *= scaleX;
	_ymin *= scaleY;
	_ymax *= scaleY;
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
	scaleAllParticlePosition(scaleX, scaleY);
}

void PeriodicControlWorld::updatePeriodicBoundary_stressControlUniform(const double& stress)
{
	double hydroStress = _totalStress.trace() / 2.0;
	double deltaVel = (stress - hydroStress) * _xlen * _dt;
	_velX += deltaVel;
	_velY = _velX;
	double scale = (_xmin + deltaVel * _dt) / _xmin;
	_xmin += deltaVel * _dt; 
	_xmax = -_xmin;
	_xlen = _xmax - _xmin;
	_ymin += deltaVel * _dt; 
	_ymax = -_ymin; 
	_ylen = _ymax - _ymin;
	scaleAllParticlePosition(scale, scale);
}

void PeriodicControlWorld::updatePeriodicBoundary_stressControlX(const double& stress)
{
	double xminPre = _xmin;
	double stressX = _totalStress.xx();
	if(abs(stressX) < stress) _velX *= (1.0 - _dt * 1.0);
	double deltaVel = (stress - stressX) * _xlen / _mass * _dt;
	_velX += deltaVel;
	_xmin += _velX * _dt;
	_xmax = -_xmin;
	_xlen = _xmax - _xmin;

	double scale = _xmin / xminPre;

	scaleAllParticlePosition(scale, 1.0);
}

void PeriodicControlWorld::updatePeriodicBoundary_stressControlY(const double& stress)
{
	double yminPre = _ymin;
	double stressY = _totalStress.yy();
	if (abs(stressY) < stress) _velY *= (1.0 - _dt * 1.0);
	double deltaVel = (stress - stressY) * _ylen / _mass * _dt;
	_velY += deltaVel;
	_ymin += _velY * _dt;
	_ymax = -_ymin;
	_ylen = _ymax - _ymin;

	double scale = _ymin / yminPre;

	scaleAllParticlePosition(1.0, scale);
}

//void PeriodicControlWorld::updatePeriodicBoundary_stressControlUniform(const double& stress)
//{
//	double xminPre = _xmin;
//	double yminPre = _ymin;
//	double hydroStress = _totalStress.trace() / 2.0;
//	double deltaVel = (stress - hydroStress) * _xlen * _dt;
//	_velX += deltaVel;
//	_velY = _velX;
//	_xmin += _velX * _dt;
//	_xmax = -_xmin;
//	_xlen = _xmax - _xmin;
//	_ymin += _velY * _dt;
//	_ymax = -_ymin;
//	_ylen = _ymax - _ymin;
//
//	double scale = _xmin / xminPre;
//
//	scaleAllParticlePosition(scale, scale);
//}

void PeriodicControlWorld::updatePeriodicBoundary_stressControlUniform_strainLimit(const double& stress)
{
	double xminPre = _xmin;
	double yminPre = _ymin;
	double hydroStress = _totalStress.trace() / 2.0;
	double deltaVel = (stress - hydroStress) * _xlen * _dt;
	_velX += deltaVel;
	if ((abs(_velX) / (_xlen / 2.0)) > _maxStrainRate) {
		_velX = abs(_velX) / _velX * (_xlen / 2.0) * _maxStrainRate;
	}
	_velY = _velX;
	
	_xmin += _velX * _dt;
	_xmax = -_xmin;
	_xlen = _xmax - _xmin;
	_ymin += _velY * _dt;
	_ymax = -_ymin;
	_ylen = _ymax - _ymin;

	double scale = _xmin / xminPre;

	scaleAllParticlePosition(scale, scale);
}




