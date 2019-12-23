#include "FrictionTermSolver.h"



FrictionTermSolver::FrictionTermSolver() :r(0.017)
{
}


FrictionTermSolver::~FrictionTermSolver()
{
}


void FrictionTermSolver::evaluateFrictionTermRHS(double gra, double * fphys, double * frhs)
{
	signed char *status = meshunion->status;
	int *K = meshunion->K;
	int *Np = meshunion->cell_p->Np;
	const int num = (*K)*(*Np);
	double g = gra;

	double *fphys_1 = fphys;
	double *fphys_2 = fphys + num;
	double *fphys_3 = fphys + 2 * num;
	double *frhs_2 = frhs + num;
	double *frhs_3 = frhs + 2 * num;

	//bool *ind;
	//requestmemory(&ind, K);
	double *s;
	requestmemory(&s, Np, K);

	for (int i = 0; i < (*K); i++)
	{
		if (status[i] == (signed char)enumSWERegion::Wet) {
			for (int j = 0; j < (*Np); j++) {
				s[i*(*Np) + j] = sqrt(pow(fphys_2[i*(*Np) + j], 2) +
					pow(fphys_3[i*(*Np) + j], 2)) /
					pow(fphys_1[i*(*Np) + j], 2);

				frhs_2[i*(*Np) + j] = frhs_2[i*(*Np) + j]
					- g * r*r*(fphys_2[i*(*Np) + j] * s[i*(*Np) + j]) /
					(pow(fphys_1[i*(*Np) + j], 1.0 / 3.0));

				frhs_3[i*(*Np) + j] = frhs_3[i*(*Np) + j]
					- g * r*r*(fphys_3[i*(*Np) + j] * s[i*(*Np) + j]) /
					(pow(fphys_1[i*(*Np) + j], 1.0 / 3.0));
				//ind[i] = true;
			}
			//else {
			//	ind[i] = false;
			//}
		}



		//for (int i = 0; i < (*K); i++)
		//{
		//	if (ind[i] == true)
		//	{
		//		for (int j = 0; j < (*Np); j++)
		//		{
		//			s[i*(*Np) + j] = (pow(fphys_2[i*(*Np) + j], 2) +
		//				pow(fphys_3[i*(*Np) + j], 2)) /
		//				pow(fphys_1[i*(*Np) + j], 2);

		//			frhs_2[i*(*Np) + j] = frhs_2[i*(*Np) + j]
		//				- g * r*r*(fphys_2[i*(*Np) + j] * s[i*(*Np) + j]) /
		//				(pow(fphys_1[i*(*Np) + j], 1 / 3));

		//			frhs_3[i*(*Np) + j] = frhs_3[i*(*Np) + j]
		//				- g * r*r*(fphys_3[i*(*Np) + j] * s[i*(*Np) + j]) /
		//				(pow(fphys_1[i*(*Np) + j], 1 / 3));
		//		}
		//	}
		//}



		//freememory(&ind);

	};
	freememory(&s);
}