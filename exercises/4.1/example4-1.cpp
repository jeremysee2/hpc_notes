#include <iostream>
#include <cstdlib>
using namespace std;
#define F77NAME(x) x##_
extern "C" {
    double F77NAME(ddot) (const int& n,
                          const double *x, const int& incx,
                          const double *y, const int& incy);
}

int main() {
    // Allocate memory for vectors
    const int N = 64;
    double* x = new double[N];
    double* y = new double[N];
    double z = 0.0;

    // Randomize x,y inputs
    std::srand(time(0));
    for (int i=0; i<N; i++) {
        x[i] = (double)(rand()) / RAND_MAX;
        y[i] = (double)(rand()) / RAND_MAX;
    }

    // Calling BLAS C API, prefixed with `cblas_`
    z = F77NAME(ddot)(N,x,1,y,1);

    std::cout << "Dot product result: " << z << endl;

    // Free dynamic memory
    delete[] x;
    delete[] y;
}