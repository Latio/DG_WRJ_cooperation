#pragma once
#include"NdgQuadFreeStrongFormAdvSolver2d.h"
#include"AbstractOutputFile.h"
#include<vector>
//#include<algorithm>
#include<fstream>
#include<iomanip>

extern "C" {

}
class NdgPhysMat :public SWEPreBlanaced2d
{
public:
	NdgPhysMat();
	~NdgPhysMat();

	void matSolver();
	void matEvaluateSSPRK22();
	void UpdateExternalField(double tloc, double *fphys);
	void EvaluateRHS(double *fphys, double *frhs);
	void UpdateOutputResult(double& time, double *fphys, int Nvar);
	void matEvaluateLimiter(double *fphys);
	//double* EvaluatePostFunc(double *fphys);
	//double UpdateTimeInterval(double *fphys);
	//void EvaluateSourceTerm(double *fphys);

protected:

	SWEAbstract2d sweabstract2d;
	SWEConventional2d sweconventional2d;
	NdgQuadFreeStrongFormAdvSolver2d ndgquadfreestrongformadvsolver2d;
	AbstractOutputFile abstractoutputfile;

	//double *fext;
	//int *outputfile;
	//double *limiter;
	//double *davectionSolver;
	//double *viscositySolver;
	//double *NonhydrostaticSolver;
	//double ftime;
	std::string casename;

	double *frhs;
	double *fext;
	double *fphys0;
	double *fphys;
	double *zGrad;
	std::vector<double> tidal;
	std::vector<int> obeindex;

	double tidalinterval;

	double ftime;
	int outputIntervalNum;
	int *Np;
	int *K;
	int Nfield;
	int *Nv;
	int Nvar;
	int *boundarydge_Nfp;
	int *boundarydge_Ne;
	//double gra;
	//double hmin;
	//int Np;// dimension 
	//int K;// dimension 
	//int Nfield;// dimension 
	//int Nvar;// dimension 

};

