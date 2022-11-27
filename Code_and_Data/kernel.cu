#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

__global__ void Plus(float A[], float B[], float C[], int n)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    C[i] = A[i] + B[i];
}

__global__ void Naive(char* pat, char* txt, int* match, int pattern_size, int text_size)
{
    for (int i = 0; i <= text_size - pattern_size; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < pattern_size; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == pattern_size) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            printf("Pattern found at index : %d\n", i);
    }
}



/*use pattern to compare starting with every possible position*/
__global__ void brute_force(char* text, char* pattern, int* match, int pattern_size, int text_size) {

    /*get the absolute pid*/
    int pid = threadIdx.x + blockIdx.x * blockDim.x;

    if (pid <= text_size - pattern_size) {

        int flag = 1;
        for (int i = 0; i < pattern_size; i++) {
            if (text[pid + i] != pattern[i]) {
                flag = 0;
            }
        }
        match[pid] = flag;
        if (flag)
            printf("Found pattern at index %d ", pid);
    }
}

/*use pattern to compare starting with every possible position*/
__global__ void brute_forces(char* text, char* pattern, bool match, int pattern_size, int text_size) {

    /*get the absolute pid*/
    int pid = threadIdx.x + blockIdx.x * blockDim.x;

    if (pid <= text_size - pattern_size) {

        int flag = 1;
        for (int i = 0; i < pattern_size; i++) {
            if (text[pid + i] != pattern[i]) {
                flag = 0;
            }
        }
        match = flag;
        if (flag)
            printf("Found pattern at index %d ", pid);
    }
}
int cap_division(int x, int y) {
    return (x + y - 1) / y;
}
int main()
{
    time_t start, end;
    time(&start);
    int size;

    std::cout << "Hello World!\n";
    string line;
    streampos sizes;
    char* memblock;
    ifstream myfile("text.txt", ios::in | ios::binary | ios::ate);
    if (myfile.is_open())
    {
        sizes = myfile.tellg();
        memblock = new char[sizes];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, sizes);
        myfile.close();

        cout << "the entire file content is in memory";
        //cout << memblock;
        char* text = memblock;
        char* pattern = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaab";
        int pattern_size = strlen(pattern);
        int* witness_array = (int*)malloc(sizeof(int) * cap_division(pattern_size, 2));
        size = strlen(text);
        char* dev_text;
        char* dev_pattern;
        int* dev_match;
        int* dev_output;
        int* dev_witness;
        int* dev_blockoutput;
        int number_of_threads = cap_division(pattern_size, 2);
        int number_of_blocks = (size + number_of_threads - 1) / number_of_threads;
        int* blockoutput = (int*)malloc(number_of_blocks * sizeof(int));
        cudaMalloc((void**)&dev_text, size * sizeof(char));
        cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
        cudaMalloc((void**)&dev_match, size * sizeof(int));
        //cudaMalloc((void **)&dev_output, sizeof(int)*size);
        cudaMalloc((void**)&dev_blockoutput, sizeof(int) * number_of_blocks);
        cudaMalloc((void**)&dev_witness, sizeof(int) * cap_division(pattern_size, 2));

        cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);

        cudaMemcpy(dev_text, text, size * sizeof(char), cudaMemcpyHostToDevice);

        cudaMemcpy(dev_witness, witness_array, cap_division(pattern_size, 2) * sizeof(int), cudaMemcpyHostToDevice);


        bool ma = 1;

        //brute_force << <(size + 1024 - 1) / 1024, 1024 >> > (dev_text, dev_pattern, dev_match, pattern_size, size);
        brute_forces << <(size + 1024 - 1) / 1024, 1024 >> > (dev_text, dev_pattern, ma, pattern_size, size);

        cudaFree(dev_text);
        cudaFree(dev_pattern);
        cudaFree(dev_match);
        time(&end);
        int timeuse = end - start;
        int seconds = (difftime(start, end));
        
        cout << "total time is " << timeuse << "ms" << endl;

        delete[] memblock;
    }
    else cout << "Unable to open file";
    
    return 0;
}