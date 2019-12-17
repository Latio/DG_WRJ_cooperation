#include "SWEPreBlanaced2d.h"



SWEPreBlanaced2d::SWEPreBlanaced2d()
{
}


SWEPreBlanaced2d::~SWEPreBlanaced2d()
{
}

void SWEPreBlanaced2d::EvaluateFlux(double *fphys, double *E, double *G)
{
	signed char *status = meshunion->status;
	sweprebalancevolumeflux2d.evaluate(hmin, gra, status, fphys, E, G);

};