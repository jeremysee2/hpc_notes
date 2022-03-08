#include <iostream>
#include <chrono>
#include <thread>
using namespace std;

#include <omp.h>

int main() {
    #pragma omp parallel
    {
        #pragma omp single
        {
            #pragma omp task
            {
                cout << "Hello world from thread "
                    << omp_get_thread_num() << endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            }
            #pragma omp task
            {
                cout << "Goodbye world from thread "
                    << omp_get_thread_num() << endl;
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            }
            cout << "MASTER has finished creating tasks." << endl;
        }
    }
    cout << "End of single region." << endl;
}