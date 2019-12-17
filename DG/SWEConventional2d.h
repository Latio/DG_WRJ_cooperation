#pragma once
#include "SWEAbstract2d.h"
extern "C" {
	void c_EvaluatePostFunc2d(double hmin_, double *fphys_, double *hc_, double *qxc_, double *qyc_, int *K_, int *Np_);
}
class SWEConventional2d :
	public SWEAbstract2d
{
public:
	SWEConventional2d();
	~SWEConventional2d();
	void EvaluatePostFunc(double *fphys);
	void UpdateWetDryState(double *fphys);
};

