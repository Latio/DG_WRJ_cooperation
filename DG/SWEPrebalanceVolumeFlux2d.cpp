#include "SWEPrebalanceVolumeFlux2d.h"



SWEPrebalanceVolumeFlux2d::SWEPrebalanceVolumeFlux2d()
{
}


SWEPrebalanceVolumeFlux2d::~SWEPrebalanceVolumeFlux2d()
{
}

void SWEPrebalanceVolumeFlux2d::evaluate(double hmin, double gra, signed char *status, double *fphys, double *E, double *G)
{
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;

	c_EvaluateFlux2d(hmin, gra, status, fphys, Np, K, E, G);
};