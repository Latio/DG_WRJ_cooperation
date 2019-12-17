#include "SWETopographySourceTerm2d.h"



SWETopographySourceTerm2d::SWETopographySourceTerm2d()
{
}


SWETopographySourceTerm2d::~SWETopographySourceTerm2d()
{
}

void SWETopographySourceTerm2d::EvaluateTopographySourceTerm(double gra, double *fphys, double *zGrad, double *frhs)
{
	signed char *status = meshunion->status;
	int *Np = meshunion->cell_p->Np;
	int *K = meshunion->K;
	int Nfield = meshunion->Nfield;
	int Nvar = 3;

	double *frhs_temp;
	requestmemory(&frhs_temp, Np, K, Nvar);
	c_EvaluateSourceTopography2d(gra, status, fphys, zGrad, frhs_temp, Np, K, Nfield);


	//for (int i = 0; i < 1080 * 3; i++)
	//{
	//	std::cout << i << " ;before " << frhs_temp[i] << std::endl;
	//}
	const int num = (*Np)*(*K)*Nvar;

	cblas_daxpy(num, 1, frhs_temp, 1, frhs, 1);//????wrong????
	//for (int i = 0; i < 1080 * 3; i++)
	//{
	//	std::cout << i << " ; " << frhs[i] << std::endl;
	//}

	freememory(&frhs_temp);
};
