#pragma once

#include <math.h>
#include "MathVector.h"

class matrix2d
{
private:
	double XX, XY, YX, YY;

public:
	matrix2d() : XX(0.0), XY(0.0), YX(0.0), YY(0.0) {}

	matrix2d(const double& xx, const double& xy, const double& yx, const double& yy) :
		XX(xx), XY(xy), YX(yx), YY(yy)
	{}

	inline double xx() { return XX; }
	inline double xy() { return XY; }
	inline double yx() { return YX; }
	inline double yy() { return YY; }
	inline const double trace() const { return (XX + YY); }
	inline const double det() const { return (XX * YY - XY * YX); }

	inline matrix2d& operator+=(const matrix2d& other);
	inline matrix2d& operator/=(const double& factor);
	inline matrix2d operator*(const double& factor) const;
	inline matrix2d& operator*=(const double& factor);

	inline vector2d dot(const vector2d& v);
	inline static matrix2d dyadic(const vector2d& v1, const vector2d& v2);
	inline vector2d eigen(matrix2d& eigenVector);
	inline const void print() const { printf("\t\t%le\t%le\n\t\t%le\t%le\n", XX, XY, YX, YY); }
};

inline matrix2d& matrix2d::operator+=(const matrix2d& other)
{
	XX += other.XX;
	XY += other.XY;
	YX += other.YX;
	YY += other.YY;
	return *this;
}

inline matrix2d& matrix2d::operator/=(const double& factor)
{
	XX /= factor;
	XY /= factor;
	YX /= factor;
	YY /= factor;
	return *this;
}

inline matrix2d matrix2d::operator*(const double& factor) const
{
	return matrix2d(XX * factor, XY * factor, YX * factor, YY * factor);
}

inline matrix2d& matrix2d::operator*=(const double& factor)
{
	XX *= factor;
	XY *= factor;
	YX *= factor;
	YY *= factor;
	return *this;
}


inline vector2d matrix2d::dot(const vector2d& v) {
	return vector2d(XX * v.x() + XY * v.y(), YX * v.x() + YY * v.y());
}

inline matrix2d matrix2d::dyadic(const vector2d& v1, const vector2d& v2)
{
	return matrix2d(v1.x() * v2.x(), v1.x() * v2.y(), v1.y() * v2.x(), v1.y() * v2.y());
}

inline vector2d matrix2d::eigen(matrix2d& eigenVector)
{
	double trace = this->trace();
	double det = this->det();
	double delta = sqrt(trace * trace - 4.0 * det);
	double eigenValue1 = (trace + delta) / 2.0;
	double eigenValue2 = (trace - delta) / 2.0;
	vector2d eigenValue;
	if (abs(eigenValue1) > abs(eigenValue2)) {
		eigenValue.setX(eigenValue1);
		eigenValue.setY(eigenValue2);
	}
	else {
		eigenValue.setX(eigenValue2);
		eigenValue.setY(eigenValue1);
	}
}


