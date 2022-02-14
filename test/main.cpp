#include <iostream>
#include <cmath>
#include <iomanip>
#include "cblas.h"

using namespace std;

/**
 * @brief Fill matrix A in symmetrical banded format.
 * 
 */
void FillMatrixA(int N, double* A, int dimA) {
    // Initialise first part
    A[1] = 1.0;

    // 1x Diagonal and 1x upper diagonal
    for (int i=1; i<N; i++) {
        A[i*dimA]       = 0.5;
        A[i*dimA + 1]   = 1.0;
    }
}

/**
 * @brief Iteratively find the solution to Ax = b using Richardson method.
 */
int main() {
    // Initialise variables
    int N           = 20;
    double alpha    = 0.9;
    int dimA        = 2;

    // Allocate memory for vectors and matrices
    double* y       = new double[N];
    double* b       = new double[N];
    double* b_ax    = new double[N];
    double* x       = new double[N];
    double* A       = new double[N*dimA];

    // Task 1: Generate a vector of length y
    for (int i=0; i<N; i++) {
        y[i] = i;
    }

    // Task 2: Generate a matrix A using symmetrical banded format
    FillMatrixA(N, A, dimA);

    // Task 3: Compute the vector b = Ay
    cblas_dsbmv(CblasColMajor, CblasUpper, N, 1, 1.0, A, dimA,
                    y, 1, 0.0, b, 1);

    // Task 4: Solve Ax = b, to compute x.
    // Initialise x to some random value
    srand(time(0));
    for (int i=0; i<N; i++) {
        x[i] = double(rand()/RAND_MAX);
    }

    // Perform x(k+1) = x(k) + alpha*(b - Ax(k)). BLAS functions only.
    for (int k=0; k<1000; k++) {
        // Let b_ax = b
        cblas_dcopy(N, b, 1, b_ax, 1);
        
        // Calculate b_ax = b_ax-Ax
        cblas_dsbmv(CblasColMajor, CblasUpper, N, 1, -1.0, A, dimA,
                        x, 1, 1.0, b_ax, 1);

        // Calculate alpha*b_ax
        cblas_dscal(N, alpha, b_ax, 1);

        // Calculate x = x + alpha*b_ax
        cblas_daxpy(N, 1.0, b_ax, 1, x, 1);
    }

    // Task 5: Calculate error, the norm of ||b - Ax|| using BLAS
    // Calculate b - Ax
    cblas_dsbmv(CblasColMajor, CblasUpper, N, 1, -1.0, A, dimA,
                    x, 1, 1.0, b, 1);

    // Calculate norm of ||b - Ax||
    double norm = cblas_dnrm2(N, b, 1);
    cout << "Error: " << norm << endl;
}