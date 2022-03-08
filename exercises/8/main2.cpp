#include <iostream>
#include <chrono>
#include <thread>
using namespace std;

#include <omp.h>

void DoWork1(int p) {
    int x = rand() % 10;
    std::this_thread::sleep_for(std::chrono::milliseconds(100*x));
    cout << "Completed work1 index " << p
         << " on thread " << omp_get_thread_num() << endl;
}

void DoWork2(int p) {
    int x = rand() % 10;
    std::this_thread::sleep_for(std::chrono::milliseconds(100*x));
    cout << "Completed work2 index " << p
         << " on thread " << omp_get_thread_num() << endl;
}

int main() {
    int n;
    std::string nstr;
    cout << "Enter number of times to run work jobs: ";
    cin >> nstr;
    n = stoi(nstr);

    int i = 0;
    #pragma omp parallel
    {
        #pragma omp single
        {
            for (i=0; i<n; ++i) {
                #pragma omp task firstprivate(i)
                DoWork1(i);
            }

            #pragma omp taskwait
            for (i=0; i<n; ++i) {
                #pragma omp task firstprivate(i)
                DoWork2(i);
            }
            cout << "Finished generating tasks" << endl;
        }
        cout << "Ended single region" << endl;
    }
}