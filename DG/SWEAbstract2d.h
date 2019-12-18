#pragma once
//#include "NdgPhysMat.h"
#include "SWEFaceFluxSolver2d.h"
#include"SWEHLLNumFluxSolver2d.h"
#include"SWEPrebalanceVolumeFlux2d.h"
#include"SWETopographySourceTerm2d.h"
#include"SWEElevationLimiter2d.h"

extern "C" {
	//void surfluxSolver_evaluate(double hmin_, double gra_, double *nx_, double *ny_, MeshUnion *mesh_, InnerEdge *edge_);//(hmin, gra, nx, ny, fm, mesh, edge)
	double c_UpdateTimeInterval2d(double hmin_, double gra_, int N_, double *dx_, signed char *status_, double * const fphys_, int *Np_, int *K_, int Nfield_);
	void c_ImposeBoundaryCondition(double gra_, double *nx_, double *ny_, double *fp_, double *fext_, signed char *ftype_, int *Nfp_, int* Ne_, int Nfield_);
	void c_HydrostaticReconstruction(double hmin_, double *fm_, double *fp_, const int *Nfp_, const int *Ne_, const int Nfield_);
}

class SWEAbstract2d
{

public:

	const double gra;
	const double hmin;
	const double cfl;
	const int Nfield;
	const double Nvar;
	double *dx;

	enum enumSWERegion {
		Sponge = 3, // % sponge cell
		Wet,		//well cell(SWE)
		Dry,		//dry cell(SWE)
		PartialWet,
		PartialWetFlood,
		PartialWetDamBreak
	} enumsweregion;



	void EvaluateSurfFlux(double *nx, double *ny, double *fm, double *fluxM, int *Nfp, int *Ne);
	void EvaluateSurfNumFlux(double *nx, double *ny, double *fm, double *fp, double *fluxS, int *Nfp, int *Ne);
	void ImposeBoundaryCondition(double *nx, double *ny, double *fm, double *fp, double *fext);
	void EvaluateSourceTerm(double *fphys, double *frhs, double *zGrad);
	double UpdateTimeInterval(double *fphys);

	SWEAbstract2d();
	~SWEAbstract2d();

	SWEFaceFluxSolver2d swefacefluxsolver2d;
	SWEHLLNumFluxSolver2d swehllnumfluxsolver2d;
	SWEPrebalanceVolumeFlux2d sweprebalancevolumeflux2d;
	SWETopographySourceTerm2d swetopographysourceterm2d;
	SWEElevationLimiter2d sweelevationlimiter2d;
};

