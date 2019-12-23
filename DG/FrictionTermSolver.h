#pragma once
#include"MeshUnion.h"
#include<math.h>

extern const MeshUnion *meshunion;
extern MeshUnion mesh;

class FrictionTermSolver
{
public:
	FrictionTermSolver();
	~FrictionTermSolver();

	void evaluateFrictionTermRHS(double gra, double *fphys, double *frhs);

	double r;
	//double g;

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet = 4,		//well cell(SWE)
		Dry = 5,		//dry cell(SWE)
		PartialWet = 6,
		PartialWetFlood = 7,
		PartialWetDamBreak = 8
	} enumsweregion;
};