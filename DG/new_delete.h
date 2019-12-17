#pragma once

void requestmemory(bool **meshpoint, int*dim1);
void requestmemory(double **meshpoint, int&dim1, int&dim2);
void requestmemory(double **meshpoint, int*dim1, int*dim2);
void requestmemory(double **meshpoint, int*dim1, int*dim2, int*dim3);
void requestmemory(double **meshpoint, int*dim1, int*dim2, int dim3);
void requestmemory(double **meshpoint, int*dim1, int dim2, int dim3);
void requestmemory(double **meshpoint, int&dim1);
void requestmemory(double **meshpoint, int*dim1);

void requestmemory(int **meshpoint, int*dim1, int*dim2);
void requestmemory(int **meshpoint, int&dim1, int&dim2);
void requestmemory(int **meshpoint, int&dim1);

void freememory(double **meshpoint);

void freememory(int **meshpoint);

void freememory(signed char **meshpoint);
void freememory(bool **meshpoint);