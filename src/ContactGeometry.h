#pragma once

#include "MathVector.h"

struct ContactGeometry
{
public:
	bool isContacted = false;
	vector2d contactPoint;
	vector2d contactNormal;
	vector2d branchVectorFromP1;
	vector2d branchVectorFromP2;
	double normalOverlap = 0.0;
	vector2d relVel;
	vector2d normalRelVel;
	vector2d shearRelVel;

	ContactGeometry() {}

	void print() {
		if (isContacted) {
			printf("Contact is detected...\n");
			printf("%+25s : (%lf, %lf)\n", "Contact point", contactPoint.x(), contactPoint.y());
			printf("%+25s : (%lf, %lf)\n", "Contact normal", contactNormal.x(), contactNormal.y());
			printf("%+25s : (%lf, %lf)\n", "Branch vector from p1", branchVectorFromP1.x(), branchVectorFromP1.y());
			printf("%+25s : (%lf, %lf)\n", "Branch vector from p1", branchVectorFromP2.x(), branchVectorFromP2.y());
			printf("%+25s : (%lf, %lf)\n", "Relative velocity", relVel.x(), relVel.y());
			printf("%+25s : (%lf, %lf)\n", "Normal relative velocity", normalRelVel.x(), normalRelVel.y());
			printf("%+25s : (%lf, %lf)\n", "Shear relative velocity", shearRelVel.x(), shearRelVel.y());
		}
	}
};