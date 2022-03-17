#include "MathMatrix2d.h"
#include "GluedParticle.h"


#define NUM_GRID_TO_EVALUATE_PARTICLE_SHAPE 100

void GluedParticle::calculateParticleProperties(const double& density)
{
	printf("[MESSAGE] Function called...GluedParticle::calculateParticleProperties\n");
	printf("          %s :\n", __FILE__);
	printf("          %s : %d\n", __FUNCTION__, __LINE__);

	double xmin = 1e10, xmax = -1e10;
	double ymin = 1e10, ymax = -1e10;
	double rmax = 0.0;
	for (auto& p : _subParticles) {
		printf("[DEBUG] sub particle ID = %d\n", p->getID());
		printf("[DEBUG] position = (%le, %le)\n", p->getPosX(), p->getPosY());
		xmin = (p->getPosX() - p->getOuterRadius()) < xmin ? (p->getPosX() - p->getOuterRadius()) : xmin;
		xmax = (p->getPosX() + p->getOuterRadius()) > xmax ? (p->getPosX() + p->getOuterRadius()) : xmax;
		ymin = (p->getPosY() - p->getOuterRadius()) < ymin ? (p->getPosY() - p->getOuterRadius()) : ymin;
		ymax = (p->getPosY() + p->getOuterRadius()) > ymax ? (p->getPosY() + p->getOuterRadius()) : ymax;
		rmax = p->getOuterRadius() > rmax ? p->getOuterRadius() : rmax;
		printf("=========================\n");
	}
	printf("[DEBUG] domain = (%le, %le), (%le, %le)\n", xmin, xmax, ymin, ymax);
	printf("[DEBUG] max radius = %le\n", rmax);
	xmin -= 0.1 * rmax; xmax += 0.1 * rmax; ymin -= 0.1 * rmax; ymax += 0.1 * rmax;
	double lx = xmax - xmin, ly = ymax - ymin;
	double dl = (lx < ly ? lx : ly) / NUM_GRID_TO_EVALUATE_PARTICLE_SHAPE;
	int nx = static_cast<int>(floor(lx / dl));
	int ny = static_cast<int>(floor(ly / dl));
	double dlx = lx / nx;
	double dly = ly / ny;
	vector<double> xx(nx), yy(ny);
	for (int i = 0; i < nx; ++i) xx[i] = xmin + dlx * (i + 0.5);
	for (int j = 0; j < ny; ++j) yy[j] = ymin + dly * (j + 0.5);
	double dvolume = dlx * dly;
	double dmass = dvolume * density;
	vector2d massCenter(0.0, 0.0);
	double volume = 0.0, mass = 0.0, inertia = 0.0;
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (auto& p : _subParticles) {
				if (p->pointInsideParticle(vector2d(xx[i], yy[j]))) {
					massCenter += vector2d(xx[i], yy[j]) * dmass;
					volume += dvolume;
					mass += dmass;
					break;
				}
			}
		}
	}
	massCenter /= mass;
	_volume = volume;
	_mass = mass;
	_invMass = 1.0 / _mass;
	_position = massCenter;
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (auto& p : _subParticles) {
				if (p->pointInsideParticle(vector2d(xx[i], yy[j]))) {
					inertia += dmass * _position.getDistance(vector2d(xx[i], yy[j]));
					break;
				}
			}
		}
	}
	_inertia = inertia;
	_invInertia = 1.0 / inertia;
	for (auto& sp : _subParticles) {
		_branchVectors.push_back(sp->getPos() - _position);
		sp->setMass(_mass);
		sp->setInertia(_inertia);
	}

	printf("          mass center = (%le, %le)\n", _position.x(), _position.y());
	printf("          mass        = %le\n", _mass);
	printf("          inertia     = %le\n", _inertia);
}

bool GluedParticle::calculateContactGeometry(
	const BaseParticle* other,
	ContactGeometry& contactGeometry,
	const vector2d& offset,
	const vector2d& offsetVel)
{
	return false;
}

void GluedParticle::takeTimeIntegral(
	const double& dt,
	vector<vector2d>::iterator iterForce,
	vector<double>::iterator iterTorque)
{

	//printf("=========================================================\n");
	//printf("[DEBUG] Function called...GluedParticle::takeTimeIntegral\n");
	//printf("=========================================================\n");

	vector2d totalForce(0.0, 0.0);
	double totalTorque = 0.0;
	for (auto& bv : _branchVectors) {
		totalForce += *(++iterForce);
		totalTorque += *(++iterTorque);
		totalTorque += bv.cross(*iterForce);
	}
	_velocity = (_velocity * (1.0 - _GlobalDamping_ * dt) + totalForce * _invMass * dt)
		/ (1.0 + _GlobalDamping_ * dt);
	_position += _velocity * dt;
	_rotationVel = (_rotationVel * (1.0 - _GlobalDamping_ * dt) + totalTorque * _invInertia * dt)
		/ (1.0 + _GlobalDamping_ * dt);
	_rotation += _rotationVel * dt;
	double angle = _rotationVel * dt;
	matrix2d R(std::cos(angle), -std::sin(angle),
			   std::sin(angle),  std::cos(angle));
	for (int i = 0; i < _subParticles.size(); ++i) {
		_branchVectors[i] = R.dot(_branchVectors[i]);
		_subParticles[i]->setVel(_velocity + _branchVectors[i].cross(_rotationVel));
		_subParticles[i]->setPos(_position + _branchVectors[i]);
		_subParticles[i]->setRotVel(_rotationVel);
		_subParticles[i]->accumulateRot(angle);
	}

}

const matrix2d GluedParticle::getStressFromContactList() const
{
	matrix2d stress(0.0, 0.0, 0.0, 0.0);
	for (int i = 0; i < _subParticles.size(); ++i) {
		auto& particle = _subParticles[i];
		auto& branchVector = _branchVectors[i];
		for (auto& iter : particle->getContactInfo()) {
			const ContactInformation& contactInfo = iter.second;
			stress += matrix2d::dyadic(contactInfo._force, (contactInfo._branchVector + branchVector));
		}
	}
	return stress;
}

const bool GluedParticle::isRattlerParticle() const
{
	for (auto& sp : _subParticles) {
		if (!sp->isRattlerParticle()) return false;
	}
	return true;
}

void GluedParticle::generateMatlabCodeForCheck() const
{
	printf("==================================================================\n");
	printf("+++++++++++++++++++++ Matlab Code for check ++++++++++++++++++++++\n");
	printf("==================================================================\n\n");
	printf("clear;clc;\n");
	printf("x = [ ");
	for (auto& sp : _subParticles) printf("%le ", sp->getPosX());
	printf("];\n");
	printf("y = [ ");
	for (auto& sp : _subParticles) printf("%le ", sp->getPosY());
	printf("];\n");
	printf("r = [ ");
	for (auto& sp : _subParticles) printf("%le ", sp->getOuterRadius());
	printf("];\n");
	printf("bx = [ ");
	for (auto& bv : _branchVectors) printf("%le ", bv.x());
	printf("];\n");
	printf("by = [ ");
	for (auto& bv : _branchVectors) printf("%le ", bv.y());
	printf("];\n");
	printf("mcx = %le;\n", _position.x());
	printf("mcy = %le;\n", _position.y());
	printf("figure;\nhold on;\naxis equal;\n");
	printf("for i = 1:length(r)\n");
	printf("\trectangle(\"Position\", [x(i)-r(i),y(i)-r(i),2*r(i),2*r(i)], \"Curvature\",1);\n");
	printf("\tquiver(mcx,mcy,bx(i),by(i),0,\"Color\",\"b\");\n");
	printf("end\n");
	printf("scatter(x,y,20,\"xr\")\n");
	printf("scatter(mcx,mcy,20,\"or\")\n");
	printf("\n==================================================================\n");
}

void GluedParticle::print2File_shape(vector<FILE*>& files)
{
	print2File_radius(files[0]);
	print2File_aspectRatio(files[1]);
}
