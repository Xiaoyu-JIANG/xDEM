#include "EllipticParticle.h"


using std::min;


void EllipticParticle::calculateParticleProperties(const double& density)
{
	_volume = M_PI * _radius * _radius / _aspectRatio;
	_mass = _volume * density;
	_invMass = 1.0 / _mass;
	_inertia = _mass * (_radius * _radius + _radius * _radius / _aspectRatio / _aspectRatio) / 4.0;
	_invInertia = 1.0 / _inertia;
}

bool EllipticParticle::calculateContactGeometry(
	const BaseParticle* otherParticle, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	if (const EllipticParticle* other = dynamic_cast<const EllipticParticle*>(otherParticle)) {
		if (!takeCoarseDetection(otherParticle, offset)) {
			contactGeometry.isContacted = false;
			return false;
		}
		return contactWithEllipticParticle(other, contactGeometry, offset, offsetVel);
	}
	else {
		printf("Todo: more particle types...\n"); exit(-1);
	}
}

bool EllipticParticle::contactWithEllipticParticle(
	const EllipticParticle* other, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	/* 
	The mean idea comes from
	https://journals.aps.org/pre/abstract/10.1103/PhysRevE.75.061709
	Distance of closest approach of two arbitrary hard ellipses in two dimensions
	Xiaoyu Zheng and Peter Palffy-Muhoray
	*/
	double temp;
	// Start accurate detection
	// Define basic properties
	double a1 = _radius, a2 = other->getRadius();
	double b1 = a1 / _aspectRatio, b2 = a2 / other->getAspectRatio();
	double z1 = _rotation, z2 = other->getRot();
	vector2d d(other->getPos()+offset-_position);
    vector2d d_hat = d / d.getLength();
    vector2d k1(cos(z1),sin(z1));
    vector2d k2(cos(z2),sin(z2));
    matrix2d k1k1 = matrix2d::dyadic(k1,k1);
    matrix2d k2k2 = matrix2d::dyadic(k2,k2);
    // Eccentricity of ellipses
    double e1 = sqrt(1.0 - pow2(b1)/pow2(a1));
    double e2 = sqrt(1.0 - pow2(b2)/pow2(a2));
    // Unit matrix
    matrix2d I(1.0, 0.0, 0.0, 1.0);
    // Transformation matrix: elps1 -> circ1: eq(4)
    matrix2d T = (I + k1k1*(b1/a1-1.0)) / b1;
    // Eta: eq(6)
    double eta = a1/b1 - 1.0;
    // Branch vector after transformation: eq(10)
    vector2d d_prime_hat = T * d_hat;
	vector2d d_prime0 = T * d;
    // Coefficient matrix of elps2 after transformation: eq(8)
    matrix2d A_prime = (I+k1k1*eta) * (I-k2k2*pow2(e2)) * (I+k1k1*eta) * pow2(b1)/pow2(b2);
    // Eigenvalues and EigenVectors of A_prime eq(9)
    vector2d k_prime_plus, k_prime_minus;
    vector2d lambda = A_prime.eigen(k_prime_plus, k_prime_minus);
    double a2_prime = 1.0 / sqrt(lambda.y());
    double b2_prime = 1.0 / sqrt(lambda.x());
	if (d_prime0.getLength() > 1.0 + a2_prime) {
		contactGeometry.isContacted = false;
		return false;
	}
	if (abs(vector2d::dot(d_prime0,k_prime_plus)) > 1.0 + b2_prime) {
		contactGeometry.isContacted = false;
		return false;
	}
    // !! make sure that both vectors 'k_prime_plus' and 'k_prime_minus' have
    // positive value when being taken dot product with vector 'd_prime_hat'. In
    // such a way eq(21) are garanteed to be positive. Thus the square root of
    // eq(23) and (24) be positive.
    if (vector2d::dot(k_prime_minus,d_prime_hat) < 0) k_prime_minus = -k_prime_minus;
    if (vector2d::dot(k_prime_plus,d_prime_hat) < 0) k_prime_plus = -k_prime_plus;
    // Derive tan(theta) and delta: eq(20-21)
    double delta = pow2(a2_prime) / pow2(b2_prime) - 1;
    double tanPhi = vector2d::dot(k_prime_minus, d_prime_hat) / vector2d::dot(k_prime_plus, d_prime_hat);
    double tanPhi2 = pow2(tanPhi);
    // Coefficients of standard quartic: eq(26)
    double A = -(1.0+tanPhi2) / pow2(b2_prime);
    double B = -2.0*(1.0+tanPhi2+delta) / b2_prime;
    double C = -tanPhi2 - pow2(1.0+delta) + (1.0+(1.0+delta)*tanPhi2)/pow2(b2_prime);
    double D = 2.0*(1.0+tanPhi2)*(1.0+delta) / b2_prime;
    double E = (1.0+tanPhi2+delta) * (1.0+delta);

	/*
    // Ferrari's method: eq(28)
    double alpha = (-3.0*pow2(B))/(8.0*pow2(A)) + C/A;
    double beta = pow3(B)/(8.0*pow3(A)) - (B*C)/(2.0*pow2(A)) + D/A;
    double gamma = (-3.0*pow4(B))/(256.0*pow4(A)) + (C*pow2(B))/(16.0*pow3(A)) - (B*D)/(4.0*pow2(A)) + E/A;
    double P = -pow2(alpha)/12.0 - gamma;
    double Q = -pow3(alpha)/108.0 + alpha*gamma/3.0 - pow2(beta)/8.0;
    // Special treatment for 'U' being a complex number. However, 'q' will eventually be a real number.
    // U = Ure + Uim * i
    temp = pow2(Q)/4.0 + pow3(P)/27.0;
    double Ure = -Q/2.0;
    double Uim = sqrt(abs(temp));
    if (temp < 0) {
        MathFunc::cpow(Ure, Uim, Ure, Uim, 1.0/3.0);
    }
    else {
        Ure = Ure + Uim;
        Uim = 0.0;
        if (Ure < 0) {
            MathFunc::cpow(Ure, Uim, Ure, Uim, 1.0/3.0);
        }
        else {
            Ure = std::pow(Ure, 1.0/3.0);
        }
    }
    // Derive Y: eq(29-30)
    double Yre = 0.0, Yim = 0.0;
    if (abs(Ure) < 1e-8 && abs(Uim) < 1e-8) {
		// Comments: The principal value of the cubic root of a negative real number is a complex number.
		double Qre = Q;
		double Qim = 0.0;
		MathFunc::cpow(Qre, Qim, Qre, Qim, 1.0 / 3.0);
		Yre = -5.0 / 6.0 * alpha - Qre;
		Yim = -Qim;
    }
    else {
        Yre = -5.0/6.0 * alpha + Ure - P/3.0 * Ure/(pow2(Ure)+pow2(Uim));
        Yim = Uim + P/3.0 * Uim/(pow2(Ure)+pow2(Uim));
    }
    // Derive q: eq(31)
    //// Define alpha+2y=q0
    double q = 0.0;
    double q0re = alpha + 2.0*Yre;
    double q0im = 2.0*Yim;
    MathFunc::cpow(q0re, q0im, q0re, q0im, 0.5);
    double q1re = - (3.0*alpha + 2.0*Yre + 2.0*beta * q0re/(pow2(q0re)+pow2(q0im)));
    double q1im = - (            2.0*Yim - 2.0*beta * q0im/(pow2(q0re)+pow2(q0im)));
    MathFunc::cpow(q1re, q1im, q1re, q1im, 0.5);
    q = -B/(4.0*A) + 0.5*(q0re + q1re);
	if (abs(alpha + 2.0 * Yre) < 1e-4 && abs(2.0 * Yim) < 1e-4) {
		q = -B / (4.0 * A) + sqrt((-alpha + sqrt(pow2(alpha) - 4.0 * gamma)) / 2.0);
	}
	*/

	/*
	Instead of the original method proposed in the paper
	Newton-Raphson Method is used to get 'q'.
	Because it was figured out that the desired root 'q' is the positive real one,
	The initial guess of 'q' is chosen as sqrt(delta+1.0), because according to eq(23) the value of 'q'
	ranges from 0 to sqrt(delta+1.0)
	*/
	double q = NewtonRaphsonForQuartic(A, B, C, D, E, sqrt(delta+1.0));


    // Closest distance of two ellipses after transformation: eq(33)
    temp = (pow2(q) - 1.0) / delta;
	double temp0 = temp;
	double d_prime = 0.0;
	if (temp < 0) temp = 0.0;
	if (delta == 0 || !isfinite(temp)) {
		d_prime = 1.0 + a2_prime;
	}
	else {
		d_prime = sqrt(temp*pow2((1.0+b2_prime*(1.0+delta)/q)) + (1.0-temp)*pow2(1.0+b2_prime/q));
	}
    // Closest distance of two ellipses before transformation: eq(36)
    double dist = d_prime / sqrt(1.0 - pow2(vector2d::dot(k1,d_hat))*pow2(e1)) * b1;
    
    // Calculate contact information
    // [0] Pre-calculate overlap
    double overlap0 = dist - d.getLength();
    if (overlap0 < 0.0) {
		contactGeometry.isContacted = false;
		return false;
	}
    // [1] other information
	if (delta == 0) {
		vector2d tangentPoint = T.inv() * d_hat;
		vector2d contactPoint = tangentPoint * d.getLength() / dist;
		matrix2d A1 = (I - k1k1 * pow2(e1)) / pow2(b1);
		vector2d contactNorm = (A1 * tangentPoint);
		contactNorm.normalize();
		vector2d overlap = d_hat * overlap0;
		overlap = contactNorm * overlap.dot(contactNorm);
		contactGeometry.isContacted = true;
		contactGeometry.contactPoint = contactPoint + _position;
		contactGeometry.contactNormal = contactNorm;
		contactGeometry.normalOverlap = overlap.getLength();
		contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - _position;
		contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (other->getPos() + offset);
	}
	else {
		double sinPsi2 = (pow2(q) - 1.0) / delta;
		if (sinPsi2 < 0) sinPsi2 = 0.0;
		double cosPsi2 = 1.0 - sinPsi2;
		if (cosPsi2 < 0) cosPsi2 = 0.0;
		matrix2d K_prime_inv(k_prime_minus.x(), k_prime_plus.x(), k_prime_minus.y(), k_prime_plus.y());
		vector2d tangentPoint = T.inv() * (K_prime_inv * vector2d(sqrt(sinPsi2), sqrt(cosPsi2)));
		vector2d contactPoint = tangentPoint * d.getLength() / dist;
		matrix2d A1 = (I - k1k1 * pow2(e1)) / pow2(b1);
		vector2d contactNorm = (A1 * tangentPoint);
		contactNorm.normalize();
		vector2d overlap = d_hat * overlap0;
		overlap = contactNorm * overlap.dot(contactNorm);
		contactGeometry.isContacted = true;
		contactGeometry.contactPoint = contactPoint + _position;
		contactGeometry.contactNormal = contactNorm;
		contactGeometry.normalOverlap = overlap.getLength();
		contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - _position;
		contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (other->getPos() + offset);

		/*if (contactGeometry.checkDataCorrectness()) {
		 	std::cout << "a1 = " << a1 << std::endl;
		 	std::cout << "b1 = " << b1 << std::endl;
		 	std::cout << "z1 = " << z1 << std::endl;
		 	std::cout << "a2 = " << a2 << std::endl;
		 	std::cout << "b2 = " << b2 << std::endl;
		 	std::cout << "z2 = " << z2 << std::endl;
		 	std::cout << "d = "; d.print();
			std::cout << "P = " << P << std::endl;
			std::cout << "Q = " << Q << std::endl;
		 	std::cout << "U = " << Ure << "\t" << Uim << std::endl;
		 	std::cout << "Y = " << Yre << "\t" << Yim << std::endl;
		 	std::cout << "q = " << q << std::endl;
		 	std::cout << "d_prime = " << d_prime << std::endl;
		 	std::cout << "dist = " << dist << std::endl;
		 	std::cout << "temp0 = " << temp0 << std::endl;
		 	std::cout << "delta = " << delta << std::endl;
		 	std::cout << "lambda = "; lambda.print();

		 	std::cout << "sinPsi2 = " << sinPsi2 << std::endl;
		 	std::cout << "cosPsi2 = " << cosPsi2 << std::endl;
		 	std::cout << "K_prime_inv = "; K_prime_inv.print();
		 	std::cout << "tangentPoint = "; tangentPoint.print();
		 	std::cout << "T = "; T.print();
		 	std::cout << "Tinv = "; T.inv().print();
		 	std::cout << "contactPoint = "; contactPoint.print();
		 	std::cout << "A1 = "; A1.print();
		 	std::cout << "contactNorm = "; contactNorm.print();
		 	std::cout << "overlap = "; overlap.print();


			double q0 = sqrt(alpha + 2.0 * Yre);
			std::cout << "q0 = " << q0 << std::endl;
			std::cout << "(-(3.0 * alpha + 2.0 * Yre + 2.0 * beta / q0)) = " 
				<< (-(3.0 * alpha + 2.0 * Yre + 2.0 * beta / q0)) << std::endl;

		 	std::cin.get();
		}*/
	}
    
	/*if (contactGeometry.checkDataCorrectness()) {
		std::cout << "a1 = " << a1 << std::endl;
		std::cout << "b1 = " << b1 << std::endl;
		std::cout << "z1 = " << z1 << std::endl;
		std::cout << "a2 = " << a2 << std::endl;
		std::cout << "b2 = " << b2 << std::endl;
		std::cout << "z2 = " << z2 << std::endl;
		std::cout << "d = "; d.print();
		std::cout << "U = " << Ure << "\t" << Uim << std::endl;
		std::cout << "Y = " << Yre << "\t" << Yim << std::endl;
		std::cout << "q = " << q << std::endl;
		std::cout << "d_prime = " << d_prime << std::endl;
		std::cout << "dist = " << dist << std::endl;
		std::cout << "temp0 = " << temp0 << std::endl;
		std::cout << "delta = " << delta << std::endl;
		std::cout << "lambda = "; lambda.print();

		std::cin.get();
	}*/
	
	calculateRelativeVelocityAtContactPoint(other, contactGeometry, offsetVel);
	return true;
}


bool EllipticParticle::contactWithEllipticParticle0(
	const EllipticParticle* other, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	// Start accurate detection
	double ai = _radius, aj = other->getRadius();
	double bi = ai / _aspectRatio, bj = aj / other->getAspectRatio();
	double xi = _position.x(), xj = other->getPos().x() + offset.x();
	double yi = _position.y(), yj = other->getPos().y() + offset.y();
	double zi = _rotation, zj = other->getRot();
	double szi = sin(zi), czi = cos(zi);
	double szj = sin(zj), czj = cos(zj);
	// coordination transformation
    // move the origin of ellipse i to (0,0) and rotate its major-semi axis to x-axis
	double z0 = zj - zi;
	double x0 = (xj - xi) * czi + (yj - yi) * szi;
	double y0 = -(xj - xi) * szi + (yj - yi) * czi;
	// coefficients of ellipse j equation after transformation
	double c0 = pow(cos(z0) / aj, 2) + pow(sin(z0) / bj, 2);
	double c1 = pow(sin(z0) / aj, 2) + pow(cos(z0) / bj, 2);
	double c2 = sin(2.0 * z0) / aj / aj - sin(2.0 * z0) / bj / bj;
	double c3 = -2.0 * x0 * c0 - y0 * c2;
	double c4 = -2.0 * y0 * c1 - x0 * c2;
	double c5 = x0 * x0 * c0 + y0 * y0 * c1 + x0 * y0 * c2 - 1.0;
	// coefficients of quartic equation
	double kk = bi * bi / (ai * ai);
	double p0 = kk * c2 * c2 + pow(c0 - c1 * kk, 2);
	double p1 = 2.0 * kk * c2 * c4 + 2.0 * (c0 - c1 * kk) * c3;
	double p2 = -kk * c2 * c2 * ai * ai + kk * c4 * c4 + 2.0 * (c0 - c1 * kk) * (c5 + c1 * kk * ai * ai) + c3 * c3;
	double p3 = 2.0 * c3 * (c5 + c1 * kk * ai * ai) - 2.0 * kk * c2 * c4 * ai * ai;
	double p4 = pow((c5 + c1 * kk * ai * ai), 2) - kk * c4 * c4 * ai * ai;
	// Solve the quartic polynomia using general formula.
	// WIKI: https://en.wikipedia.org/wiki/Quartic_function
	double alphaInit, betaInit, gammaInit, deltaInit;
	bool rootFoundGeneral = solveQuarticGeneral(p0, p1, p2, p3, p4,
		alphaInit, betaInit, gammaInit, deltaInit);
	if (!rootFoundGeneral) {
		contactGeometry.isContacted = false;
		return false;
	}
	p1 /= p0;
	p2 /= p0;
	p3 /= p0;
	p4 /= p0;
	p0 /= p0;
	double alpha, beta, gamma, delta;
	double finalError = fastQuarticSolver(p1, p2, p3, p4, alphaInit, betaInit, gammaInit, deltaInit,
		alpha, beta, gamma, delta);
	// solve roots by {alpha,beta}, {gamma,delta}
	std::vector<double> xroots;
	std::vector<double> yroots;
	double delta1 = alpha * alpha - 4.0 * beta;
	double delta2 = gamma * gamma - 4.0 * delta;
	if (delta1 >= delta2) {
		xroots.push_back((-alpha - sqrt(delta1)) / 2.0);
		xroots.push_back((-alpha + sqrt(delta1)) / 2.0);
	}
	else {
		xroots.push_back((-gamma + sqrt(delta2)) / 2.0);
		xroots.push_back((-gamma - sqrt(delta2)) / 2.0);
	}
	if (xroots.size() != 2) {
		contactGeometry.isContacted = false;
		return false;
	}
	if (abs(xroots[0]) / ai > 1.0 && abs(xroots[0]) / ai <= 1.0 + 1e-6) xroots[0] = copysign(ai, xroots[0]);
	if (abs(xroots[1]) / ai > 1.0 && abs(xroots[1]) / ai <= 1.0 + 1e-6) xroots[1] = copysign(ai, xroots[1]);
	// calculate corresponding y-values from x-roots
	for (std::size_t i = 0; i < xroots.size(); i++) {
		yroots.push_back(bi * sqrt(1.0 - pow(xroots[i] / ai, 2)));
		double error1 = c0 * xroots[i] * xroots[i] + c1 * yroots[i] * yroots[i] + c2 * xroots[i] * yroots[i]
			+ c3 * xroots[i] + c4 * yroots[i] + c5;
		double error2 = c0 * xroots[i] * xroots[i] + c1 * yroots[i] * yroots[i] - c2 * xroots[i] * yroots[i]
			+ c3 * xroots[i] - c4 * yroots[i] + c5;
		yroots[i] = (abs(error1) <= abs(error2)) ? yroots[i] : -yroots[i];
		if (min(abs(error1), abs(error2)) > 1e-8) {
			contactGeometry.isContacted = false;
			return false;
		}
		if (!isfinite(yroots[i])) {
			contactGeometry.isContacted = false;
			return false;
		}
		if (abs(yroots[i]) > bi) {
			contactGeometry.isContacted = false;
			return false;
		}
	}
	// find two intersections of two ellipses
	vector2d contactPoint1;
	vector2d contactPoint2;
	contactPoint1.setX(xroots[0] * czi - yroots[0] * szi + xi);
	contactPoint1.setY(xroots[0] * szi + yroots[0] * czi + yi);
	contactPoint2.setX(xroots[1] * czi - yroots[1] * szi + xi);
	contactPoint2.setY(xroots[1] * szi + yroots[1] * czi + yi);
	contactGeometry.contactPoint = (contactPoint1 + contactPoint2) / 2.0; // find real contact point
	if (!isfinite(contactGeometry.contactPoint.getLength())) {
		contactGeometry.isContacted = false;
		return false;
	}
	// find contact tangent to find contact normal
	vector2d tangent = contactPoint1 - contactPoint2;
	if (tangent.getLength() == 0.0) {
		contactGeometry.isContacted = false;
		return false;
	}
	tangent.normalize();
	contactGeometry.contactNormal = vector2d(tangent.y(), -tangent.x());
	if (!isfinite(contactGeometry.contactNormal.getLength())) {
		contactGeometry.isContacted = false;
		FILE* logFile;
		logFile = fopen("logFile.dat", "a");
		fprintf(logFile, "[Log] contact normal...\n");
		fprintf(logFile, "tangent = (%le, %le)\n", tangent.x(), tangent.y());
		fprintf(logFile, "normal  = (%le, %le)\n", contactGeometry.contactNormal.x(), contactGeometry.contactNormal.y());
		fprintf(logFile, "============================================\n\n");
		fflush(logFile);
		fclose(logFile);
		return false;
	}
	vector2d bv = contactGeometry.contactPoint - vector2d(xi, yi);
	if (contactGeometry.contactNormal.dot(bv) < 0.0) contactGeometry.contactNormal = -contactGeometry.contactNormal;
	/* [2] Find contact overlap dn */
	double k = contactGeometry.contactNormal.y() / contactGeometry.contactNormal.x();
	if (abs(k) > 100000) k = copysign(100000, k);
	if (!isfinite(k)) k = (contactGeometry.contactNormal.y() > 0) ? 100000 : -100000;
	// ellipse i
	double aai = pow(czi / ai, 2) + pow(szi / bi, 2);
	double bbi = pow(szi / ai, 2) + pow(czi / bi, 2);
	double cci = czi * szi * (1.0 / ai / ai - 1.0 / bi / bi);
	double alphai = contactGeometry.contactPoint.y() - yi - k * contactGeometry.contactPoint.x();
	double p_i1 = aai + bbi * k * k + 2.0 * cci * k;
	double p_i2 = -2.0 * aai * xi + 2.0 * bbi * k * alphai + 2.0 * cci * (alphai - k * xi);
	double p_i3 = aai * xi * xi + bbi * alphai * alphai - 2.0 * cci * alphai * xi - 1.0;
	double delta_i = sqrt(p_i2 * p_i2 - 4.0 * p_i1 * p_i3);
	double xoi1 = (-p_i2 - delta_i) / (2.0 * p_i1);
	double xoi2 = (-p_i2 + delta_i) / (2.0 * p_i1);
	double yoi1 = k * (xoi1 - contactGeometry.contactPoint.x()) + contactGeometry.contactPoint.y();
	double yoi2 = k * (xoi2 - contactGeometry.contactPoint.x()) + contactGeometry.contactPoint.y();
	double disti1 = pow(xoi1 - contactGeometry.contactPoint.x(), 2) + pow(yoi1 - contactGeometry.contactPoint.y(), 2);
	double disti2 = pow(xoi2 - contactGeometry.contactPoint.x(), 2) + pow(yoi2 - contactGeometry.contactPoint.y(), 2);
	double disti = (disti1 > disti2) ? sqrt(disti2) : sqrt(disti1);
	// ellipse j
	double aaj = pow(czj / aj, 2) + pow(szj / bj, 2);
	double bbj = pow(szj / aj, 2) + pow(czj / bj, 2);
	double ccj = czj * szj * (1.0 / aj / aj - 1.0 / bj / bj);
	double alphaj = contactGeometry.contactPoint.y() - yj - k * contactGeometry.contactPoint.x();
	double p_j1 = aaj + bbj * k * k + 2.0 * ccj * k;
	double p_j2 = -2.0 * aaj * xj + 2.0 * bbj * k * alphaj + 2.0 * ccj * (alphaj - k * xj);
	double p_j3 = aaj * xj * xj + bbj * alphaj * alphaj - 2.0 * ccj * alphaj * xj - 1.0;
	double delta_j = sqrt(p_j2 * p_j2 - 4.0 * p_j1 * p_j3);
	double xoj1 = (-p_j2 - delta_j) / (2.0 * p_j1);
	double xoj2 = (-p_j2 + delta_j) / (2.0 * p_j1);
	double yoj1 = k * (xoj1 - contactGeometry.contactPoint.x()) + contactGeometry.contactPoint.y();
	double yoj2 = k * (xoj2 - contactGeometry.contactPoint.x()) + contactGeometry.contactPoint.y();
	double distj1 = pow(xoj1 - contactGeometry.contactPoint.x(), 2) + pow(yoj1 - contactGeometry.contactPoint.y(), 2);
	double distj2 = pow(xoj2 - contactGeometry.contactPoint.x(), 2) + pow(yoj2 - contactGeometry.contactPoint.y(), 2);
	double distj = (distj1 > distj2) ? sqrt(distj2) : sqrt(distj1);
	contactGeometry.normalOverlap = disti + distj;
	// A modification for large k
	if (!isfinite(disti)) contactGeometry.normalOverlap = 2.0 * distj;
	if (!isfinite(distj)) contactGeometry.normalOverlap = 2.0 * disti;
	contactGeometry.isContacted = true;
	contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - _position;
	contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (other->getPos() + offset);
	calculateRelativeVelocityAtContactPoint(other, contactGeometry, offsetVel);
	if (!isfinite(contactGeometry.normalOverlap)) {
		contactGeometry.isContacted = false;
		FILE* logFile;
		logFile = fopen("logFile.dat", "a");
		fprintf(logFile, "[Log] contact normal overlap...\n");
		fprintf(logFile, "k = %le\n", k);
		fprintf(logFile, "tangent = (%le, %le)\n", tangent.x(), tangent.y());
		fprintf(logFile, "normal  = (%le, %le)\n", contactGeometry.contactNormal.x(), contactGeometry.contactNormal.y());
		fprintf(logFile, "aai = %le\tbbi = %le\tcci = %le\n", aai, bbi, cci);
		fprintf(logFile, "alphai = %le\n", alphai);
		fprintf(logFile, "p_i1 = %le\tp_i2 = %le\tp_i3 = %le\n", p_i1, p_i2, p_i3);
		fprintf(logFile, "delta_i = %le\n", delta_i);
		fprintf(logFile, "(xoi1, yoi1) = (%le, %le)\n", xoi1,yoi1);
		fprintf(logFile, "(xoi2, yoi2) = (%le, %le)\n", xoi2,yoi2);
		fprintf(logFile, "disti   = %le\t%le\t%le\n", disti1, disti2, disti);
		fprintf(logFile, "aaj = %le\tbbj = %le\tccj = %le\n", aaj, bbj, ccj);
		fprintf(logFile, "alphaj = %le\n", alphaj);
		fprintf(logFile, "p_j1 = %le\tp_j2 = %le\tp_j3 = %le\n", p_j1, p_j2, p_j3);
		fprintf(logFile, "delta_j = %le\n", delta_j);
		fprintf(logFile, "(xoj1, yoj1) = (%le, %le)\n", xoj1, yoj1);
		fprintf(logFile, "(xoj2, yoj2) = (%le, %le)\n", xoj2, yoj2);
		fprintf(logFile, "distj   = %le\t%le\t%le\n", distj1, distj2, distj);
		fprintf(logFile, "============================================\n\n");
		fflush(logFile);
		fclose(logFile);
		return false;
	}
	return true;
}

const bool EllipticParticle::pointInsideParticle(const vector2d& point) const
{
	double cosa = cos(_rotation);
	double sina = sin(_rotation);
	double dd = (_radius / _aspectRatio) * (_radius / _aspectRatio);
	double DD = _radius * _radius;
	double a = pow(cosa * (point.x() - _position.x()) + sina * (point.y() - _position.y()), 2);
	double b = pow(sina * (point.x() - _position.x()) - cosa * (point.y() - _position.y()), 2);
	return (a / dd + b / DD) <= 1.0;
}

void EllipticParticle::print2File_shape(vector<FILE*>& files)
{
	print2File_radius(files[0]);
	print2File_aspectRatio(files[1]);
}

#define FORMAT_EllipticParticle_print2File_particle \
"%7d\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%4d\n"
void EllipticParticle::print2File_particle(FILE* file)
{
	fprintf(file, FORMAT_EllipticParticle_print2File_particle,
		_id, _outerRadius, _radius, _aspectRatio,
		_position.x(), _position.y(),
		_velocity.x(), _velocity.y(),
		_rotation, _rotationVel,
		_volume, _mass, _invMass, _inertia, _invInertia,
		static_cast<int>(_contactInfo.size()));
	for (auto& iter : _contactInfo) {
		iter.second.print2File_contactInfo(file);
	}
}

// WIKI: https://en.wikipedia.org/wiki/Quartic_function
// ax ^ 4 + bx ^ 3 + cx ^ 2 + dx + e = 0 (with real coefficients and 'a' is not equal to 0)
bool EllipticParticle::solveQuarticGeneral(
	double& a, double& b, double& c, double& d, double& e, // input parameters
	double& alphaInit, double& betaInit, double& gammaInit, double& deltaInit) // output parameters
{
	// (1) calculate its discriminant delta
	// if the discriminant(delta) is negative,
	// then the equation has two distinct real rootsand two complex conugate non - real roots
	// Our object is to find two intersection points, coresponding to two real roots only
	bool rootFound = false;
	double delta0 = c * c - 3.0 * b * d + 12.0 * a * e;
	double delta1 = 2.0 * c * c * c - 9.0 * b * c * d + 27.0 * b * b * e + 27.0 * a * d * d - 72.0 * a * c * e;
	double delta = (delta1 * delta1 - 4.0 * delta0 * delta0 * delta0) / (-27.0);
	if (delta <= 1e-12 && delta > 0.0) {
		delta = 0.0;
	}
	else if (delta > 1e-12) {
		return rootFound;
	}
	// MAIN IDEA !!
	// if any negative value in any square root, stop the calculation
	double p = (8.0 * a * c - 3.0 * b * b) / (8.0 * a * a);
	double q = (b * b * b - 4.0 * a * b * c + 8.0 * a * a * d) / (8.0 * a * a * a);
	double t = (delta1 + sqrt(delta * (-27.0))) / 2.0;
	t = pow(t, 1. / 3.);
	double temp = (-2.0 / 3.0) * p + (t + delta0 / t) / (3.0 * a);
	if (temp <= 0.0) return rootFound;
	double s = 0.5 * sqrt(temp);
	alphaInit = -(-b / (2.0 * a) - 2.0 * s);
	betaInit = (-b / (4.0 * a) - s) * (-b / (4.0 * a) - s) - 0.25 * (-4.0 * s * s - 2.0 * p + q / s);
	gammaInit = -(-b / (2.0 * a) + 2.0 * s);
	deltaInit = (-b / (4.0 * a) + s) * (-b / (4.0 * a) + s) - 0.25 * (-4.0 * s * s - 2.0 * p - q / s);
	rootFound = true;
	return rootFound;
}

// The fast quartic solver
	// Peter Strobach(2010)
	// https://www.sciencedirect.com/science/article/pii/S0377042710002128
double EllipticParticle::fastQuarticSolver(
	double& a, double& b, double& c, double& d,
	double& alphaInit, double& betaInit, double& gammaInit, double& deltaInit,
	double& alpha, double& beta, double& gamma, double& delta)
{
	// Initialization of set {alpha,beta,gamma,delta}
	alpha = alphaInit;
	beta = betaInit;
	gamma = gammaInit;
	delta = deltaInit;
	// Initial error from initialized {alpha,beta,gamma,delta} to real ones
	double e1 = a - alpha - gamma;
	double e2 = b - beta - alpha * gamma - delta;
	double e3 = c - beta * gamma - alpha * delta;
	double e4 = d - beta * delta;
	double e_old = abs(e1) + abs(e2) + abs(e3) + abs(e4);
	double e_new = 1.0;
	// Iteration start
	int counter = 0;
	double u23, u33, L43, u44;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	while (abs((e_old - e_new) / e_old) >= 1e-3) {
		e_old = e_new;
		counter++;
		if (counter > 20) break;
		u23 = alpha - gamma;
		u33 = beta - delta - gamma * u23;
		L43 = -delta * u23 / u33;
		u44 = beta - delta - L43 * u23;;
		x1 = e1;
		x2 = e2 - gamma * x1;
		x3 = e3 - delta * x1 - gamma * x2;
		x4 = e4 - delta * x2 - L43 * x3;
		y4 = x4 / u44;
		y3 = (x3 - u23 * y4) / u33;
		y2 = x2 - u23 * y3 - y4;
		y1 = x1 - y3;
		alpha = alpha + y1;
		beta = beta + y2;
		gamma = gamma + y3;
		delta = delta + y4;
		e1 = a - alpha - gamma;
		e2 = b - beta - alpha * gamma - delta;
		e3 = c - beta * gamma - alpha * delta;
		e4 = d - beta * delta;
		e_new = abs(e1) + abs(e2) + abs(e3) + abs(e4);
		if (e_new < 1e-8) break;
	}
	return e_new;
}



double EllipticParticle::NewtonRaphsonForQuartic(const double A, const double B, const double C, const double D, const double E, const double init)
{
	double A0 = 4.0 * A;
	double B0 = 3.0 * B;
	double C0 = 2.0 * C;
	double r0 = init;
	double r = 0.0;
	double tol = 1e-8;
	while (abs(r0 - r) > tol) {
		r = r0;
		r0 -= (A * pow4(r0) + B * pow3(r0) + C * pow2(r0) + D * r0 + E) / (A0 * pow3(r0) + B0 * pow2(r0) + C0 * r0 + D);
	}
	return r0;
}