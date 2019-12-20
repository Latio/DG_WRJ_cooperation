#include "NdgPhysMat.h"
#define max(a, b) ((a > b) ? a : b)

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

NdgPhysMat::NdgPhysMat() :frhs(NULL),
startTime(0),
finalTime(259200),
outputIntervalNum(500),
tidalinterval(600)/*潮流数据间隔*/,
abstractoutputfile("..//..//20191208.nc", 259200.0 / 500.0, 500)
{
	Np = meshunion->cell_p->Np;
	K = meshunion->K;
	boundarydge_Nfp = meshunion->boundarydge_p->Nfp;
	boundarydge_Ne = meshunion->boundarydge_p->Ne;
	//Nfield = meshunion->Nfield;
	//Nvar = 3;

	requestmemory(&fphys, Np, K, Nfield);
	requestmemory(&fext, boundarydge_Nfp, boundarydge_Ne, Nfield);
	requestmemory(&zGrad, Np, K, 2);

	netCDF::NcFile dataFile("init_fphys.nc", netCDF::NcFile::read);
	netCDF::NcVar fphys_v = dataFile.getVar("fphys");
	fphys_v.getVar(fphys);
	netCDF::NcVar zGrad_v = dataFile.getVar("zGrad");
	zGrad_v.getVar(zGrad);

	double *ind, *temp_ftoe1, *bot;
	requestmemory(&ind, boundarydge_Nfp, boundarydge_Ne);
	requestmemory(&temp_ftoe1, boundarydge_Ne, boundarydge_Nfp);
	requestmemory(&bot, Np, K);

	/////////////////////////////////////////////////
	for (int i = 0; i < *boundarydge_Nfp; i++)
	{
		cblas_dcopy(*boundarydge_Ne, meshunion->boundarydge_p->FToE, 2, temp_ftoe1 + i, *boundarydge_Nfp);
	}
	cblas_daxpy((*boundarydge_Ne)*(*boundarydge_Nfp), *Np, temp_ftoe1, 1, ind, 1);
	cblas_daxpy((*boundarydge_Ne)*(*boundarydge_Nfp), 1, meshunion->boundarydge_p->FToN1, 1, ind, 1);
	for (int i = 0; i < (*boundarydge_Nfp)*(*boundarydge_Ne); i++)
	{
		ind[i] = ind[i] - (*Np) - 1;
	}

	double *fext_4 = fext + 3 * (*boundarydge_Ne)*(*boundarydge_Nfp);
	double *fphys_4 = fphys + 3 * (*Np)*(*K);
	cblas_dcopy((*Np)*(*K), fphys_4, 1, bot, 1);
	for (int i = 0; i < (*boundarydge_Nfp)*(*boundarydge_Ne); i++)
	{
		fext_4[i] = bot[(int)ind[i]];
	}
	freememory(&ind);
	freememory(&temp_ftoe1);
	freememory(&bot);
	////////////////////////////////////////////////

	typedef enum {
		NdgEdgeInner = 0,
		NdgEdgeGaussEdge = 1,
		NdgEdgeSlipWall = 2,
		NdgEdgeNonSlipWall = 3,
		NdgEdgeZeroGrad = 4,
		NdgEdgeClamped = 5,
		NdgEdgeClampedDepth = 6,
		NdgEdgeClampedVel = 7,
		NdgEdgeFlather = 8,
		NdgEdgeNonLinearFlather = 9,
		NdgEdgeNonLinearFlatherFlow = 10,
		NdgEdgeNonReflectingFlux = 11
	} NdgEdgeType;

	signed char *ftype = meshunion->boundarydge_p->ftype;
	for (int i = 0; i < *boundarydge_Ne; i++)
	{
		NdgEdgeType type = (NdgEdgeType)ftype[i];
		if (ftype[i] == NdgEdgeClampedDepth)
		{
			obeindex.push_back(i);
		}
	}

	ifstream data("TideElevation.txt");//read tidal data
	if (!data.is_open())
	{
		std::cout << "Error File Path !!!" << std::endl;
		system("pause");
	}
	double point_tidal;
	while (data >> point_tidal)
		tidal.push_back(point_tidal);


	data.close();

}


NdgPhysMat::~NdgPhysMat()
{
	freememory(&fphys);
	freememory(&fext);
	freememory(&zGrad);

	std::cout << "析构NdgPhyMat" << std::endl;
}


void NdgPhysMat::matSolver()
{
	matEvaluateSSPRK22();
}


void NdgPhysMat::matEvaluateSSPRK22()
{
	clock_t begintime, endtime;
	begintime = clock();

	const int num = (*K)*(*Np)*Nvar;

	double time = startTime;
	double ftime = finalTime;

	double *fphys0;
	requestmemory(&fphys0, Np, K, Nvar);

	abstractoutputfile.ncFile_create(Np, K, Nvar);

	while (time < ftime)
	{
		double dt = UpdateTimeInterval(fphys)*0.4;
		cout << dt << endl;
		if (time + dt > ftime) {
			dt = ftime - time;
		}

		cblas_dcopy(num, fphys, 1, fphys0, 1);

		for (int intRK = 0; intRK < 2; intRK++) {

			double tloc = time + dt;
			UpdateExternalField(tloc, fphys);

			requestmemory(&frhs, Np, K, Nvar);
			EvaluateRHS(fphys, frhs);

			cblas_daxpy(num, dt, frhs, 1, fphys, 1);


			EvaluateLimiter(fphys);

			std::ofstream out("D:\\Desktop\\input.txt");
			if (!out)
			{
				std::cerr << "open error!" << std::endl;
			}

			for (int j = 0; j < 433; j++) {

				out << j + 1;
				for (size_t i = 0; i < 100; i++)
				{
					out << "    " << tidal[j * 100 + i] << "\n";
				}
				out << "\n";

			}
			cout << "************************************************\n";

			EvaluatePostFunc(fphys);//Update status



			freememory(&frhs);
		}

		cblas_dscal(num, 0.5, fphys, 1);
		cblas_daxpy(num, 0.5, fphys0, 1, fphys, 1);

		time = time + dt;
		UpdateOutputResult(time, fphys, Nvar);


		double timeRatio = time / ftime;
		std::cout << "____________________finished____________________: " << timeRatio << std::endl;
	}

	endtime = clock();
	std::cout << "\n\nRunning Time : " << endtime - begintime << " ms\n" << endl;

	freememory(&fphys0);
}


void NdgPhysMat::EvaluateRHS(double *fphys, double *frhs)
{
	ndgquadfreestrongformadvsolver2d.evaluateAdvectionRHS(fphys, frhs, fext);
	sweabstract2d.EvaluateSourceTerm(fphys, frhs, zGrad);
};

//void NdgPhysMat::UpdateOutputResult(double time, double *fphys) {};
void NdgPhysMat::UpdateExternalField(double tloc, double *fphys)
{
	const int benfp = *meshunion->boundarydge_p->Nfp;
	const int bene = *meshunion->boundarydge_p->Ne;
	const int obnum = benfp * obeindex.size();

	const double delta = tidalinterval;

	const int s1 = ceil(tloc / delta);//double s1 = floor(tloc / delta) + 1;
	//const int s2 = s1 + 1;
	const double alpha1 = (delta*s1 - tloc) / delta;
	double alpha2 = (tloc - delta * (s1 - 1)) / delta;

	std::vector<double> fnT;

	for (int i = 0; i < obnum; i++) {
		double temp = tidal[(s1 - 1)*obnum + i] * alpha1 + tidal[s1*obnum + i] * alpha2;
		fnT.push_back(temp);
	}

	double *fext_4 = fext + 3 * benfp * bene;
	for (int i = 0; i < obeindex.size(); i++) {
		for (int j = 0; j < benfp; j++)
		{
			fext[obeindex[i] * benfp + j] = max(fnT[i*benfp + j] - fext_4[obeindex[i] * benfp + j], 0);
		}
	}
}


void NdgPhysMat::UpdateOutputResult(double& time, double *fphys, int Nvar)
{
	abstractoutputfile.outputIntervalResult(time, fphys, Nvar, Np, K);
};

void NdgPhysMat::EvaluateLimiter(double *fphys)
{
	sweabstract2d.sweelevationlimiter2d.apply(fphys);
};
