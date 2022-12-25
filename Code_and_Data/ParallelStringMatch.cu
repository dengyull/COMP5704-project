
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <emmintrin.h>
#include <smmintrin.h>
#include <random>
using namespace std;
#define XSIZE 4200
#define SIGMA 256
#include <sys/time.h>

ofstream outfile;
unsigned long long todval (struct timeval *tp) {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}
__global__ void NaiveSearch(char* pat, char* txt, int* match, int pattern_size, int text_size)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    if (pid <= text_size - pattern_size) {
        //printf("pid %d ", pid);
        //int loop = 0;
        int flag = 1;
        for (int i = 0; i < pattern_size; i++) {
            if (txt[pid + i] != pat[i]) {
                flag = 0;
            }
            //loop++;
        }
        if (flag) {
            //printf("Found pattern at index %d ", pid);
            match[pid] = 1;
        }
        //printf("loop %d ", loop);
    }
}

int dobrute(char* text, char* pattern) {

    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* matches = new int[sizes];
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms1 = (todval(&t2)-todval(&t1))/1000;
    cout << "pre time: " << runtime_ms1/1000.0 << endl;
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    //cout << (sizes + 1024 - 1) / 1024 << endl;
    gettimeofday(&t1,0);
    NaiveSearch << <(sizes + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes);
    cudaDeviceSynchronize();
    gettimeofday(&t2,0);
    unsigned long long runtime_ms2 = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms2/1000.0 << endl;
    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);
    unsigned long long runtime_ms3 = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms3/1000.0 << endl;

    outfile << runtime_ms1 << ";" << runtime_ms2 << ";" << runtime_ms3 << ";";
    cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}

__global__ void SliceKMP(char* pat, char* txt, int* match, int patlen, int txtlen, int* lps, int range)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    int i = 0 + pid * range;
    int j = 0;
    while (((txtlen - i) >= (patlen - j)) && (i < range + pid * range + patlen -1)){
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == patlen) {
            match[(i - j)] = 1;
            j = lps[j - 1];
        }
        else if (i < txtlen && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
}

__global__ void KMPSearch(char* pat, char* txt, int* match, int patlen, int txtlen, int* lps)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    printf("pid %d ", pid);
    int i = 0;
    int j = 0;

    while ((txtlen - i) >= (patlen - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == patlen) {
            match[(i - j)] = 1;
            j = lps[j - 1];
        }
        else if (i < txtlen && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
}
void LPSArray(char* pat, int M, int* lps)
{
    int len = 0;
    lps[0] = 0;
    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else
        {
            if (len != 0) {
                len = lps[len - 1];
            }
            else
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void doKMP(char* text, char* pattern) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* dev_lps;
    int* matches = new int[sizes];
    for (int i = 0; i < sizes; i++) {
        //matches[i] = 0;
    }
    int* lps = new int[pattern_size];
    LPSArray(pattern, pattern_size, lps);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_lps, pattern_size * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lps, lps, pattern_size * sizeof(int), cudaMemcpyHostToDevice);
    cout << pattern_size << endl;
    cout << sizes << endl;
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    KMPSearch << <1, 1 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, dev_lps);
    cudaDeviceSynchronize();

    gettimeofday(&t2,0);
    unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms/1000.0 << endl;
    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms/1000.0;

    cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
}


int dosliceKMP(char* text, char* pattern, int range) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* dev_lps;
    int* matches = new int[sizes];
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    int* lps = new int[pattern_size];
    LPSArray(pattern, pattern_size, lps);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_lps, pattern_size * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lps, lps, pattern_size * sizeof(int), cudaMemcpyHostToDevice);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms1 = (todval(&t2)-todval(&t1))/1000;
    cout << "pre time: " << runtime_ms1/1000.0 << endl;
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    int tid = sizes / range;
    //cout << tid << endl;
    gettimeofday(&t1,0);
    SliceKMP << <(tid + 1024 - 1) /1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, dev_lps, range);
    
    cudaDeviceSynchronize();
    gettimeofday(&t2,0);
    unsigned long long runtime_ms2 = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms2/1000.0 << endl;
    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);
    unsigned long long runtime_ms3 = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms3/1000.0 << endl;
    outfile << runtime_ms1 << ";" << runtime_ms2 << ";" << runtime_ms3 << ";";

    //cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}


void suffixes(char* x, int m, int* suff) {
    int f, g, i;
    suff[m - 1] = m;
    g = m - 1;
    for (i = m - 2; i >= 0; --i) {
        if (i > g && suff[i + m - 1 - f] < i - g)
            suff[i] = suff[i + m - 1 - f];
        else {
            if (i < g)
                g = i;
            f = i;
            while (g >= 0 && x[g] == x[g + m - 1 - f])
                --g;
            suff[i] = f - g;
        }
    }
}
void preBmGs(char* x, int m, int bmGs[]) {
    int i, j, suff[XSIZE];
    suffixes(x, m, suff);
    for (i = 0; i < m; ++i) bmGs[i] = m;
    j = 0;
    for (i = m - 1; i >= 0; --i)
        if (suff[i] == i + 1)
            for (; j < m - 1 - i; ++j)
                if (bmGs[j] == m)
                    bmGs[j] = m - 1 - i;
    for (i = 0; i <= m - 2; ++i)
        bmGs[m - 1 - suff[i]] = m - 1 - i;
}

void preBmBc(char* x, int m, int bmBc[]) {
    int i;
    for (i = 0; i < SIGMA; ++i)
        bmBc[i] = m;
    for (i = 0; i < m - 1; ++i)
        bmBc[x[i]] = m - i - 1;
}
__global__ void bm(char* x, char* y, int* match, int patlen, int txtlen, int bmGs[], int bmBc[])
{
    int i, j;
    /* Searching */
    j = 0;
    while (j <= txtlen - patlen) {
        for (i = patlen - 1; i >= 0 && x[i] == y[i + j]; --i);
        if (i < 0) {
            match[j] = 1;
            j += bmGs[0];
        }
        else
            j += max(bmGs[i], bmBc[y[i + j]] - patlen + 1 + i);
    }
}
__global__ void slicebm(char* x, char* y, int* match, int patlen, int txtlen, int bmGs[], int bmBc[], int range)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    int i, j;
    /* Searching */
    j = 0 + pid * range;
    while ((j <= txtlen - patlen) && (j < range + pid * range + patlen - 1)) {
        for (i = patlen - 1; i >= 0 && x[i] == y[i + j]; --i);
        if (i < 0) {
            match[j] = 1;
            j += bmGs[0];
        }
        else
            j += max(bmGs[i], bmBc[y[i + j]] - patlen + 1 + i);
    }
}

void dobm(char* text, char* pattern) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* dev_lps;
    int* matches = new int[sizes];
    for (int i = 0; i < sizes; i++) {
        //matches[i] = 0;
    }
    int bmGs[XSIZE], bmBc[SIGMA];
    int* dev_bmGs;
    int* dev_bmBc;
    preBmGs(pattern, pattern_size, bmGs);
    preBmBc(pattern, pattern_size, bmBc);
    int* lps = new int[pattern_size];
    LPSArray(pattern, pattern_size, lps);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_bmGs, XSIZE * sizeof(int));
    cudaMalloc((void**)&dev_bmBc, SIGMA * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_bmGs, bmGs, XSIZE * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_bmBc, bmBc, SIGMA * sizeof(int), cudaMemcpyHostToDevice);
    cout << pattern_size << endl;
    cout << sizes << endl;
    bm << <1, 1 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, dev_bmGs, dev_bmBc);

    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }

    cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
}
int doslicebm(char* text, char* pattern, int range) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* dev_lps;
    int* matches = new int[sizes];
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    int bmGs[XSIZE], bmBc[SIGMA];
    int* dev_bmGs;
    int* dev_bmBc;
    preBmGs(pattern, pattern_size, bmGs);
    preBmBc(pattern, pattern_size, bmBc);
    int* lps = new int[pattern_size];
    LPSArray(pattern, pattern_size, lps);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_bmGs, XSIZE * sizeof(int));
    cudaMalloc((void**)&dev_bmBc, SIGMA * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_bmGs, bmGs, XSIZE * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_bmBc, bmBc, SIGMA * sizeof(int), cudaMemcpyHostToDevice);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms1 = (todval(&t2)-todval(&t1))/1000;
    cout << "pre time: " << runtime_ms1/1000.0 << endl;
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    int tid = sizes / range;
    //cout << tid << endl;
    
    gettimeofday(&t1,0);
    slicebm << <(tid + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, dev_bmGs, dev_bmBc, range);

    cudaDeviceSynchronize();
    gettimeofday(&t2,0);
    unsigned long long runtime_ms2 = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms2/1000.0 << endl;
    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);

    unsigned long long runtime_ms3 = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms3/1000.0 << endl;
    outfile << runtime_ms1 << ";" << runtime_ms2 << ";" << runtime_ms3 << ";";
    //cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}


int preSo(char* x, int m, unsigned int S[]) {
    unsigned int j, lim;
    int i;
    for (i = 0; i < 256; ++i)
        S[i] = ~0;
    for (lim = i = 0, j = 1; i < m; ++i, j <<= 1) {
        S[x[i]] &= ~j;
        lim |= j;
    }
    lim = ~(lim >> 1);
    return(lim);
}


__global__ void sol(char* x, char* y, int* match, int patlen, int txtlen, unsigned int lim, unsigned int S[256])
{
    unsigned int D, k, h, p_len;
    int j, count;
    p_len = patlen;
    patlen = 32;
    /* Searching */
    for (D = ~0, j = 0; j < txtlen; ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            k = 0;
            h = j - patlen + 1;
            while (k < p_len && x[k] == y[h + k]) k++;
            if (k == p_len) {
                match[j - patlen + 1] = 1;
            }
        }
    }
}
__global__ void slicesol(char* x, char* y, int* match, int patlen, int txtlen, unsigned int lim, unsigned int S[256], int range)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int D, k, h, p_len;
    int j, count;
    p_len = patlen;
    patlen = 32;
    /* Searching */
    for (D = ~0, j = 0 + pid * range; (j < txtlen && (j < range + pid * range + patlen - 1)); ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            k = 0;
            h = j - patlen + 1;
            while (k < p_len && x[k] == y[h + k]) k++;
            if (k == p_len) {
                match[j - patlen + 1] = 1;
            }
        }
    }
}

__global__ void so(char* x, char* y, int* match, int patlen, int txtlen, unsigned int lim, unsigned int S[256])
{
    unsigned int D;
    int j;
    for (D = ~0, j = 0; j < txtlen; ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            match[j - patlen + 1] = 1;
        }
    }
}
__global__ void sliceso(char* x, char* y, int* match, int patlen, int txtlen, unsigned int lim, unsigned int S[256],int range)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int D;
    int j;
    for (D = ~0, j = 0 + pid * range; (j < txtlen && (j < range + pid * range + patlen - 1)); ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            match[j - patlen + 1] = 1;
        }
    }
}
int doso(char* text, char* pattern) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* matches = new int[sizes];
    for (int i = 0; i < sizes; i++) {
        //matches[i] = 0;
    }
    unsigned int lim;
    unsigned int S[256];
    unsigned int* dev_S;
    lim = preSo(pattern, pattern_size, S);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_S, 256 * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_S, S, 256 * sizeof(int), cudaMemcpyHostToDevice);
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    if (pattern_size > 32) {
        sol << <1, 1 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, lim, dev_S);

    }
    else {
        so << <1, 1 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, lim, dev_S);

    }

    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }

    //cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}
int dosliceso(char* text, char* pattern, int range) {
    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* matches = new int[sizes];
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    unsigned int lim;
    unsigned int S[256];
    unsigned int* dev_S;
    lim = preSo(pattern, pattern_size, S);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMalloc((void**)&dev_S, 256 * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_S, S, 256 * sizeof(int), cudaMemcpyHostToDevice);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms1 = (todval(&t2)-todval(&t1))/1000;
    cout << "pre time: " << runtime_ms1/1000.0 << endl;
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    int tid = sizes / range;
    //cout << tid << endl;
    gettimeofday(&t1,0);
    if (pattern_size > 32) {
        slicesol << <(tid + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, lim, dev_S, range);

    }
    else {
        sliceso << <(tid + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, lim, dev_S, range);

    }
    cudaDeviceSynchronize();
    gettimeofday(&t2,0);
    unsigned long long runtime_ms2 = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms2/1000.0 << endl;

    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);
    unsigned long long runtime_ms3 = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms3/1000.0 << endl;
    outfile << runtime_ms1 << ";" << runtime_ms2 << ";" << runtime_ms3 << ";";

    //cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}


__global__ void RabinKarp2(char* x, char* y, int* match, int pattern_size, int text_size) {
    int da, hx, hy, i, j, count;
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    count = 0;
    for (da = i = 1; i < pattern_size; ++i)
        da = (da << 1);

    for (hy = hx = i = 0; i < pattern_size; ++i) {
        hx = ((hx << 1) + x[i]);
        hy = ((hy << 1) + y[i]);
    }

    j = 0;
    while (j <= text_size - pattern_size) {
        if (hx == hy) {
            int flag = 1;
            for (int i = 0; i < pattern_size; i++) {
                if (y[pid + j + i] != x[i]) {
                    flag = 0;
                }
            }
            if (flag) {
                match[pid + j] = 1;
                //printf("Found pattern at index %d ", pid + j);
            }
        }
        hy = ((((hy)-(y[pid + j]) * da) << 1) + (y[pid + j + pattern_size]));
        ++j;
    }
}


__global__ void RabinKarp(char* x, char* y, int* match, int pattern_size, int text_size, int range) {
    int da, hx, hy, i, j, count;
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    count = 0;
    for (da = i = 1; i < pattern_size; ++i)
        da = (da << 1);
    for (hy = hx = i = 0; i < pattern_size; ++i) {
        hx = ((hx << 1) + x[i]);
        hy = ((hy << 1) + y[pid * range + i]);
    }
    j = 0 + pid * range;
    while ((j <= text_size - pattern_size) && (j < range + pid * range + pattern_size - 1)) {

        //printf("pid * range %d ", pid * range);
        if (hx == hy) {
            int flag = 1;
            //printf("pattern_size %d ", pattern_size);
            for (int i = 0; i < pattern_size; i++) {
                if (y[j + i] != x[i]) {
                    flag = 0;
                    //printf("error at index %d ", pid * range);
                    //printf("error at index %d ", pid * range + j + i);
                }
            }
            if (flag) {
                //printf("Found pattern at index %d ", pid * range);
                match[j] = 1;
            }
        }
        hy = ((((hy)-(y[j]) * da) << 1) + (y[j + pattern_size]));
        ++j;
    }
}


int doRabinKarp(char* text, char* pattern) {

    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* matches = new int[sizes];
    for (int i = 0; i < sizes; i++) {
        //matches[i] = 0;
    }
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);
    //cout << pattern_size << endl;
    //cout << sizes << endl;
    //cout << (sizes + 1024 - 1) / 1024 << endl;
    RabinKarp2 << <1, 1 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes);

    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }

    //cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}

int dosliceRabinKarp(char* text, char* pattern, int range) {

    int pattern_size = strlen(pattern);
    int sizes = strlen(text);
    char* dev_text;
    char* dev_pattern;
    int* dev_match;
    int* matches = new int[sizes];
    struct timeval t1, t2;
    gettimeofday(&t1,0);
    cudaMalloc((void**)&dev_text, sizes * sizeof(char));
    cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
    cudaMalloc((void**)&dev_match, sizes * sizeof(int));
    cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_text, text, sizes * sizeof(char), cudaMemcpyHostToDevice);

    gettimeofday(&t2,0);
    unsigned long long runtime_ms1 = (todval(&t2)-todval(&t1))/1000;
    cout << "pre time: " << runtime_ms1/1000.0 << endl;
    int tid = sizes / range;
    gettimeofday(&t1,0);
    RabinKarp << <(tid + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, dev_match, pattern_size, sizes, range);

    cudaDeviceSynchronize();
    gettimeofday(&t2,0);
    unsigned long long runtime_ms2 = (todval(&t2)-todval(&t1))/1000;
    cout << "core time: " << runtime_ms2/1000.0 << endl;
    cudaMemcpy(matches, dev_match, sizes * sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t1,0);
    int count = 0;
    //cout << "count : " << count << endl;
    for (int i = 0; i < sizes; i++) {
        //cout << matches[i] << endl;
        if (matches[i] == 1) {
            count++;
        }
    }
    gettimeofday(&t2,0);
    unsigned long long runtime_ms3 = (todval(&t2)-todval(&t1))/1000;
    cout << "end time: " << runtime_ms3/1000.0 << endl;
    outfile << runtime_ms1 << ";" << runtime_ms2 << ";" << runtime_ms3 << ";";

    cout << "count : " << count << endl;
    cudaFree(dev_text);
    cudaFree(dev_pattern);
    cudaFree(dev_match);
    return count;
}




void Loop(char* pat, char* txt,int patlen, int txtlen)
{
    int count;

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= txtlen - patlen; i++) {
        if (txt[i] == pat[0]){
            count++;
        }
    }
    cout << "Loop " << count << endl;
}


void does(int n, char* text, int patsize, char* outf)
{

    //int n = 1;
    //if (argc > 1)
    //    n = atol(argv[1]);
    std::cout << "Hello World!\n";
    string line;
    streampos size;
    char* memblock;
    ifstream myfile(text, ios::in | ios::binary | ios::ate);
    if (myfile.is_open()){

        size = myfile.tellg();
        memblock = new char[size];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, size);
        myfile.close();
        cout << "the entire file content is in memory"<< endl;
    }
    //ifstream myfile3(argv[3], ios::in | ios::binary | ios::ate);

    char* pat = (char*)"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaab";
    /*new char[patsize+1];
    memcpy(pat, memblock + 1000, patsize);
    pat[patsize] = '\0';*/
    /*if (myfile3.is_open()){

        size = myfile3.tellg();
        pat = new char[size];
        myfile3.seekg(0, ios::beg);
        myfile3.read(pat, size);
        myfile3.close();
        cout << pat << endl;
    }*/
    cout << "pat gen"<< endl;
    outfile << pat << endl;

    struct timeval t1, t2;

    int counts = 0;
    int patlen = strlen(pat);
    int txtlen = strlen(memblock);
    gettimeofday(&t1,0);
    Loop(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
    cout << "Loop " << endl;
    printf("%f\n", runtime_ms/1000.0);


    gettimeofday(&t1,0);
    counts = dobrute(memblock, pat);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << "Naive;" << runtime_ms/1000.0 << ";count:" << counts << endl;
    gettimeofday(&t1,0);
    counts = doslicebm(memblock, pat,1000);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << "bm;" << runtime_ms/1000.0 << ";count:" << counts << endl;


    gettimeofday(&t1,0);
    counts = dosliceKMP(memblock, pat,1000);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << "KMPSearch;" << runtime_ms/1000.0 << ";count:" << counts << endl;


    gettimeofday(&t1,0);
    counts = dosliceso(memblock, pat,1000);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << "Shift Or;" << runtime_ms/1000.0 << ";count:" << counts << endl;


    gettimeofday(&t1,0);
    counts = dosliceRabinKarp(memblock, pat,1000);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << "RabinKarp;" << runtime_ms/1000.0 << ";count:" << counts << endl;
    delete[] memblock;
    //delete[] pat;
}

int main() {


    outfile.open("outputtext2.txt");
    //does(1000,(char*)"text2.txt",2, (char*)"otext2p2.txt");
    //does(1000,(char*)"text2.txt",4, (char*)"otext2p4.txt");
    //does(1000,(char*)"text2.txt",8, (char*)"otext2p8.txt");
    //does(1000,(char*)"text2.txt",16, (char*)"otext2p16.txt");
    //does(1000,(char*)"text2.txt",32, (char*)"otext2p32.txt");
    does(1000,(char*)"text2.txt",64, (char*)"otext2p64.txt");
}