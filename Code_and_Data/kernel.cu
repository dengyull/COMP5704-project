
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sys/time.h>
using namespace std;



unsigned long long todval (struct timeval *tp) {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}


__global__ void ShiftOrMatch(char* pSrc, int nSrcSize, char* pSubSrc, int nSubSrcSize, int tr)
{
    int i = threadIdx.x;
    int end;
    if (i >= tr - 1)
        end = nSrcSize - nSubSrcSize;
    else if (tr==1)
        end = nSrcSize - nSubSrcSize;
    else
        end = (i * (nSrcSize +1 / tr));
    long skip[256];
    memset(skip, -1, sizeof(skip));
    for (int i = 0; i < nSubSrcSize; i++)
    {
        skip[pSubSrc[i]] ^= (0x01 << i);
    }

    long mask = ~(0x01 << (nSubSrcSize - 1));
    long ds = -1;
    int nPos = (threadIdx.x * (nSrcSize / tr));
    while (nPos <= end)
    {
        ds = (ds << 0x01) | skip[pSrc[nPos]];
        if (~(ds | mask))
        {
            break;
        }
        nPos++;
    }
    //printf("Found pattern at index %d ", nPos - (nSubSrcSize - 1));
    //return nPos - (nSubSrcSize - 1);

}


__global__ void bm(char* x, char* y, int m, int n, int tr) {
    int pid = threadIdx.x + blockIdx.x * blockDim.x;
    int i, j, bmGs[4200], bmBc[256], count;
    int end;
    if (i >= tr - 1)
        end = n - m;
    else if (tr==1)
        end = n - m;
    else
        end = (threadIdx.x * (n +1 / tr));

    /* Preprocessing */
    int suff[4200];
    int f, g;
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
    for (i = 0; i < m; ++i) bmGs[i] = m;
    j = 0;
    for (i = m - 1; i >= 0; --i)
        if (suff[i] == i + 1)
            for (; j < m - 1 - i; ++j)
                if (bmGs[j] == m)
                    bmGs[j] = m - 1 - i;
    for (i = 0; i <= m - 2; ++i)
        bmGs[m - 1 - suff[i]] = m - 1 - i;
    for (i = 0; i < 256; ++i)
        bmBc[i] = m;
    for (i = 0; i < m - 1; ++i)
        bmBc[x[i]] = m - i - 1;
    

        /* Searching */
    j = (threadIdx.x * (n / tr));
    count = 0;
    while (j <= end) {
        for (i = m - 1; i >= 0 && x[i] == y[i + j]; --i);
        if (i < 0) {
            printf("Found pattern at index %d ", pid);
            j += bmGs[0];
        }
        else
            j += max(bmGs[i], bmBc[y[i + j]] - m + 1 + i);
    }
    //return count;
}





void slice(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    cout << N << endl;

    /* A loop to slide pat[] one by one */
    for (int i = 0; i < tread; i++) {

        char* part = new char[(1) * (N / tread) + 1];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread));
        part[(1) * (N / tread)] = '\0';
        //cout << part << endl;
        //cout << i * (N / tread) << endl;
        //cout << (i + 1) * (N / tread) << endl;
        //ShiftOrMatch(pat, part);
    }
}

__global__ void NaiveSearch(char* pat, char* txt, bool match, int pattern_size, int text_size)
{
    int pid = threadIdx.x + blockIdx.x * blockDim.x;

    if (pid <= text_size - pattern_size) {

        int flag = 1;
        for (int i = 0; i < pattern_size; i++) {
            if (txt[pid + i] != pat[i]) {
                flag = 0;
            }
        }
        if (flag) {
            printf("Found pattern at index %d ", pid);
            match = 1;
        }
    }
}

__global__ void RabinKarp2(char* x, char* y, int pattern_size, int text_size) {
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
        if (hx == hy ) {
            int flag = 1;
            for (int i = 0; i < pattern_size; i++) {
                if (y[pid + i] != x[i]) {
                    flag = 0;
                }
            }
            if (flag) {
                printf("Found pattern at index %d ", pid);
                count++;
            }
        }
        hy = ((((hy)-(y[j]) * da) << 1) + (y[j + pattern_size]));
        ++j;
    }
}

int main()
{

    time_t start, end;
    time(&start);
    string line;
    streampos size;
    char* memblock;
    ifstream myfile("gen10.txt", ios::in | ios::binary | ios::ate);
    if (myfile.is_open())
    {
        size = myfile.tellg();
        memblock = new char[size];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, size);
        myfile.close();

        cout << "the entire file content is in memory";
        //cout << memblock;
        char* text = memblock;
        char* pattern = "Then were the king's scribes called at that time in the third month, that is, the month Sivan, on the three and twentieth day thereof; and it was written according to all that Mordecai commanded unto the Jews, and to the lieutenants, and the deputies and rulers of the provinces which are from India unto Ethiopia, an hundred twenty and seven provinces, unto every province according to the writing thereof, and unto every people after their language, and to the Jews according to their writing, and according to their language.";
        struct timeval t1, t2;
        gettimeofday(&t1,0);
        int pattern_size = strlen(pattern);
        int sizes = strlen(text);
        char* dev_text;
        char* dev_pattern;
        int* dev_match;
        int* dev_output;
        int* dev_witness;
        int* dev_blockoutput;
        cudaMalloc((void**)&dev_text, sizes * sizeof(char));
        cudaMalloc((void**)&dev_pattern, pattern_size * sizeof(char));
        cudaMalloc((void**)&dev_match, sizes * sizeof(int));
        //cudaMalloc((void **)&dev_output, sizeof(int)*size);
        ofstream outfile;
        outfile.open("cudaout.txt");
        
        cudaMemcpy(dev_pattern, pattern, pattern_size * sizeof(char), cudaMemcpyHostToDevice);

        cudaMemcpy(dev_text, text, size * sizeof(char), cudaMemcpyHostToDevice);

        bool matches = false;

        NaiveSearch << <(sizes + 1024 - 1) / 1024, 1024 >> > (dev_pattern, dev_text, matches, pattern_size, size);

        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        printf("Navie: \n%f\n", runtime_ms/1000.0);
        outfile << endl << "Naive: " << runtime_ms/1000.0 << endl ;
        if (matches)
            printf("Found pattern ");
        else
            printf("nothing");
        //cudaFree(matching);

        gettimeofday(&t1,0);
        int M = strlen(pattern);
        int N = strlen(text);
        int ne[1024],me[1024];
        for (int i = 0; i < 1023; i++)
        {
            ne[i] = (i * (N / 1024)); 
            me[i] = (1 * (N / 1024)) + M;
        }
        ne[1023] = (1023 * (N / 1024));
        me[1023] = N - (1023 * (N / 1024));
        int* nes[1024];
        int* mes[1024];
        cudaMalloc((void**)&mes, 1024 * sizeof(int));
        cudaMalloc((void**)&nes, 1024 * sizeof(int));
        cudaMalloc((void**)&nes[0], 1024 * sizeof(int));
        cudaMemcpy(nes[0], ne, 1024 * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(nes[0], ne, 1024 * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(mes[0], me, 1024 * sizeof(int), cudaMemcpyHostToDevice);
        ShiftOrMatch << <1, 512 >> > (dev_text, size, dev_pattern, pattern_size,512);

        outfile << endl << "ShiftOrMatch: " << runtime_ms/1000.0 << endl ;
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        printf("\n%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        bm << <1, 1 >> > (dev_text, dev_pattern, size, pattern_size,512);

        outfile << endl << "bm: " << runtime_ms/1000.0 << endl ;
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        printf("\n%f\n", runtime_ms/1000.0);
        delete[] memblock;
    }
    else cout << "Unable to open file";
    return 0;
}
