
#include <iostream>

#include "src/MathVector.h"
#include "src/MathMatrix2d.h"

using namespace std;

int main()
{
    double a1 = 3.0;
    double b1 = 1.1;
    double z1 = 30.0/180.0*M_PI;
    double a2 = 2.0;
    double b2 = 0.9;
    double z2 = -20.0/180.0*M_PI;

    vector2d d(3.9, 0.1);
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
    // Coefficient matrix of elps2 after transformation: eq(8)
    matrix2d A_prime = (I+k1k1*eta) * (I-k2k2*pow2(e2)) * (I+k1k1*eta) * pow2(b1)/pow2(b2);
    // Eigenvalues and EigenVectors of A_prime eq(9)
    vector2d k_prime_plus, k_prime_minus;
    vector2d lambda = A_prime.eigen(k_prime_plus, k_prime_minus);
    double a2_prime = 1.0 / sqrt(lambda.y());
    double b2_prime = 1.0 / sqrt(lambda.x());
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
    // Ferrari's method: eq(28)
    double alpha = (-3.0*pow2(B))/(8.0*pow2(A)) + C/A;
    double beta = pow3(B)/(8.0*pow3(A)) - (B*C)/(2.0*pow2(A)) + D/A;
    double gamma = (-3.0*pow4(B))/(256.0*pow4(A)) + (C*pow2(B))/(16.0*pow3(A)) - (B*D)/(4.0*pow2(A)) + E/A;
    double P = -pow2(alpha)/12.0 - gamma;
    double Q = -pow3(alpha)/108.0 + alpha*gamma/3.0 - pow2(beta)/8.0;
    // Special treatment for 'U' being a complex number. However, 'q' will eventually be a real number.
    // U = Ure + Uim * i
    double temp = pow2(Q)/4.0 + pow3(P)/27.0;
    double Ure = -Q/2.0;
    double Uim = sqrt(abs(temp));
    bool isReal = false;
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
            isReal = true;
        }
    }
    // Derive Y: eq(29-30)
    double Yre = 0.0, Yim = 0.0;
    if (abs(Ure) < 1e-8 && abs(Uim) < 1e-8) {
        Yre = -5.0/6.0 * alpha - std::pow(Q, 1.0/3.0);
    }
    else {
        Yre = -5.0/6.0 * alpha + Ure - P/3.0 * Ure/(pow2(Ure)+pow2(Uim));
        Yim = Uim + P/3.0 * Uim/(pow2(Ure)+pow2(Uim));
    }
    // Derive q: eq(31)
    // Define alpha+2y=q0
    double q = 0.0;
    if (isReal) {
        double q0 = sqrt(alpha + 2.0*Yre);
        q = -B/(4.0*A) + 0.5*(q0 + sqrt(-(3.0*alpha + 2.0*Yre + 2.0*beta/q0)));
    }
    else {
        double q0re = alpha + 2.0*Yre;
        double q0im = 2.0*Yim;
        MathFunc::cpow(q0re, q0im, q0re, q0im, 0.5);
        double q1re = - (3.0*alpha + 2.0*Yre + 2.0*beta * q0re/(pow2(q0re)+pow2(q0im)));
        double q1im = - (            2.0*Yim - 2.0*beta * q0im/(pow2(q0re)+pow2(q0im)));
        MathFunc::cpow(q1re, q1im, q1re, q1im, 0.5);
        q = -B/(4.0*A) + 0.5*(q0re + q1re);
    }
    // Closest distance of two ellipses after transformation: eq(33)
    temp = (pow2(q) - 1.0) / delta;
    double d_prime = sqrt(temp*pow2((1.0+b2_prime*(1.0+delta)/q)) + (1.0-temp)*pow2(1.0+b2_prime/q));
    // Closest distance of two ellipses before transformation: eq(36)
    double dist = d_prime / sqrt(1.0 - pow2(vector2d::dot(k1,d_hat))*pow2(e1)) * b1;
    
    // Calculate contact information
    // [0] Pre-calculate overlap
    double overlap0 = dist - d.getLength();
    if (overlap0 < 0.0) cout << "no contact" << endl;
    // [1] contact point
    double sinPsi2 = (pow2(q) - 1.0) / delta;
    double cosPsi2 = 1.0 - sinPsi2;
	matrix2d K_prime_inv(k_prime_minus.x(),k_prime_plus.x(),k_prime_minus.y(),k_prime_plus.y());
    vector2d tangentPoint = T.inv() * (K_prime_inv * vector2d(sqrt(sinPsi2),sqrt(cosPsi2)));
    vector2d contactPoint = tangentPoint * d.getLength() / dist;
    matrix2d A1 = (I - k1k1*pow2(e1)) / pow2(b1);
    vector2d contactNorm = (A1 * tangentPoint);
    vector2d overlap = d_hat * overlap0;
    overlap = contactNorm * overlap.dot(contactNorm);
    cout << overlap << endl;
    

    return 0;
}