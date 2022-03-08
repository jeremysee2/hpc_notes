#include <iostream>
#include <cmath>
#include <thread>
using namespace std;

#include <omp.h>

double ComputePi(int it) {
    cout << "Running " << it << " iterations on thread " << omp_get_thread_num() << endl;
    // Do some work
    double estimate = 0.0;
    int n = 0;
    while (n < it) {
        estimate += pow(-1.0,n)/(2.0*n + 1.0);
        n++;
    }
    return 4.0*estimate;
}

void Print(double &a, double &b, double &c) {
    cout << "a = "<< a << ", b = " << b << ", c = " << c << endl;
}

int main() {
    double a = 0.0, b = 0.0, c = 0.0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task depend(out:a)
            a = ComputePi(1);
            #pragma omp task depend(out:b)
            b = ComputePi(10);
            #pragma omp task depend(out:c)
            c = ComputePi(10000);
            #pragma omp task depend(in:a,b,c)
            Print(a,b,c);
        }
    }
}