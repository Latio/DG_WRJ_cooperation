
/*共有方法创建和释放内存*/
#include<iostream>

void requestmemory(bool **meshpoint, int*dim1)
{
	*meshpoint = new bool[*dim1]();
};


void requestmemory(double **meshpoint, int&dim1, int&dim2)
{
	*meshpoint = new double[dim1*dim2]();
};


void requestmemory(double **meshpoint, int*dim1, int*dim2)
{
	*meshpoint = new double[(*dim1)*(*dim2)]();
};


void requestmemory(double **meshpoint, int*dim1, int*dim2, int*dim3)
{
	*meshpoint = new double[(*dim1)*(*dim2)*(*dim3)]();
};

void requestmemory(double **meshpoint, int*dim1, int dim2, int dim3)
{
	*meshpoint = new double[(*dim1)*(dim2)*(dim3)]();
};

void requestmemory(double **meshpoint, int*dim1, int*dim2, int dim3)
{
	*meshpoint = new double[(*dim1)*(*dim2)*dim3]();
};


void requestmemory(double **meshpoint, int&dim1)
{
	*meshpoint = new double[dim1]();
};

void requestmemory(double **meshpoint, int*dim1)
{
	*meshpoint = new double[*dim1]();
};


void requestmemory(int **meshpoint, int*dim1, int*dim2)
{
	*meshpoint = new int[(*dim1)*(*dim2)]();
};

void requestmemory(int **meshpoint, int&dim1, int&dim2)
{
	*meshpoint = new int[dim1*dim2];
};
void requestmemory(int **meshpoint, int&dim1)
{
	*meshpoint = new int[dim1]();
};



void freememory(double **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[](*meshpoint);
		*meshpoint = NULL;

	}

};

void freememory(int **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[](*meshpoint);
		*meshpoint = NULL;

	}

};

void freememory(signed char **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[] * meshpoint;
		*meshpoint = NULL;

	}

};

void freememory(bool **meshpoint)
{
	if (*meshpoint != NULL)
	{
		delete[] * meshpoint;
		*meshpoint = NULL;

	}

};