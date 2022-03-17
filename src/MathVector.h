#pragma once

#include <iostream>
#include "MathFunc.h"

class vector2d
{
private:
	double X = 0.0;
	double Y = 0.0;

public:
	vector2d() {}
	vector2d(const double& x, const double& y) : X(x), Y(y) {}
	vector2d(const vector2d& v) : X(v.X), Y(v.Y) {}
	~vector2d() {}

	inline const double x() const { return X; }
	inline const double y() const { return Y; }

	inline vector2d operator+(const vector2d& other) const;
	inline vector2d& operator+=(const vector2d& other);
	friend inline vector2d operator-(const vector2d& v);
	inline vector2d operator-(const vector2d& other) const;
	inline vector2d& operator-=(const vector2d& other);
	inline vector2d operator*(const double& factor) const;
	inline vector2d& operator*=(const double& factor);
	inline vector2d operator/(const double& factor) const;
	inline vector2d& operator/=(const double& factor);
	friend inline std::ostream& operator<<(std::ostream& os, const vector2d& v);
	friend inline std::istream& operator>>(std::istream& is, vector2d& v);

	inline double getLength() const;
	inline double getLengthSquared() const;
	inline double getDistance(const vector2d& other) const;
	inline double getDistanceSquared(const vector2d& other) const;
	inline static double getLength(const vector2d& v);
	inline static double getLengthSquared(const vector2d& v);
	inline static double getDistance(const vector2d& v1, const vector2d& v2);
	inline static double getDistanceSquared(const vector2d& v1, const vector2d& v2);

	inline void set(const double& x, const double& y) { X=x; Y=y; }
	inline void setX(const double& x) { X = x; }
	inline void setY(const double& y) { Y = y; }
	inline void moveX(const double& x) { X += x; }
	inline void moveY(const double& y) { Y += y; }
	inline void scaleX(const double& scale) { X *= scale; }
	inline void scaleY(const double& scale) { Y *= scale; }

	inline vector2d getUnitVector() const;
	inline void normalize();

	inline double dot(const vector2d& other) const;
	inline double cross(const vector2d& other) const;
	inline vector2d cross(const double& omega) const;
	inline static double dot(const vector2d& v1, const vector2d& v2);
	inline static double cross(const vector2d& v1, const vector2d& v2);
	inline static vector2d cross(const vector2d& v, const double& omega);

	const void print() const { printf("(%lf, %lf)\n", X, Y); }
};

inline vector2d vector2d::operator+(const vector2d& other) const
{
	return vector2d(X + other.X, Y + other.Y);
}

inline vector2d& vector2d::operator+=(const vector2d& other)
{
	X += other.X;
	Y += other.Y;
	return *this;
}

inline vector2d operator-(const vector2d& v)
{
	return vector2d(-v.X, -v.Y);
}

inline std::ostream& operator<<(std::ostream& os, const vector2d& v)
{
	os << v.X << "\t" << v.Y << "\t";
	return os;
}

inline std::istream& operator>>(std::istream& is, vector2d& v)
{
	is >> v.X >> v.Y;
	return is;
}

inline vector2d vector2d::operator-(const vector2d& other) const
{
	return vector2d(X - other.X, Y - other.Y);
}

inline vector2d& vector2d::operator-=(const vector2d& other)
{
	X -= other.X;
	Y -= other.Y;
	return *this;
}

inline vector2d vector2d::operator*(const double& factor) const
{
	return vector2d(X * factor, Y * factor);
}

inline vector2d& vector2d::operator*=(const double& factor)
{
	X *= factor;
	Y *= factor;
	return *this;
}

inline vector2d vector2d::operator/(const double& factor) const
{
	return vector2d(X / factor, Y / factor);
}

inline vector2d& vector2d::operator/=(const double& factor)
{
	X /= factor;
	Y /= factor;
	return *this;
}

inline double vector2d::getLength() const
{
	return sqrt(getLengthSquared());
}

inline double vector2d::getLengthSquared() const
{
	return X * X + Y * Y;
}

inline double vector2d::getDistance(const vector2d& other) const
{
	return getDistanceSquared(other);
}

inline double vector2d::getDistanceSquared(const vector2d& other) const
{
	return getLengthSquared(*this - other);
}

inline double vector2d::getLength(const vector2d& v)
{
	return sqrt(getLengthSquared(v));
}

inline double vector2d::getLengthSquared(const vector2d& v)
{
	return v.X * v.X + v.Y * v.Y;
}

inline double vector2d::getDistance(const vector2d& v1, const vector2d& v2)
{
	return sqrt(getDistanceSquared(v1, v2));
}

inline double vector2d::getDistanceSquared(const vector2d& v1, const vector2d& v2)
{
	return getLengthSquared(v1 - v2);
}

inline vector2d vector2d::getUnitVector() const
{
	double lengthSquared = getLengthSquared();
	if (lengthSquared != 0.0) {
		return *this / sqrt(lengthSquared);
	}
	else {
		return vector2d(0.0, 0.0);
	}
}

inline void vector2d::normalize()
{
	double lengthSquared = getLengthSquared();
	if (lengthSquared != 0.0) {
		*this /= sqrt(lengthSquared);
	}
	else {
		printf("__ERROR__normalize a zero-length vector...\n");
		exit(-1);
	}
}

inline double vector2d::dot(const vector2d& other) const
{
	return X * other.X + Y * other.Y;
}

inline double vector2d::cross(const vector2d& other) const
{
	return X * other.Y - Y * other.X;
}

inline vector2d vector2d::cross(const double& omega) const
{
	return vector2d(-omega * Y, omega * X);
}

inline double vector2d::dot(const vector2d& v1, const vector2d& v2)
{
	return v1.X * v2.X + v1.Y * v2.Y;
}

inline double vector2d::cross(const vector2d& v1, const vector2d& v2)
{
	return v1.X * v2.Y - v1.Y * v2.X;
}

inline vector2d vector2d::cross(const vector2d& v, const double& omega)
{
	return vector2d(-omega * v.Y, omega * v.X);
}
