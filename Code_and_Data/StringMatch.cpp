#include <iostream>
#include <string.h>
#include <fstream>

#define d 256 
#define XSIZE 4200
#define SIGMA 256

using namespace std;

#include <sys/time.h>
#include <cilk/cilk.h>
unsigned long long todval (struct timeval *tp) {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}

void NavieCal(char* pat, char* txt,int i,int M,int N){
    int j;
        /* For current index i, check for pattern match */
        for (j = 0; j < M; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == M) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            cout << "Pattern found at index " << i << endl;
}

void PNaive(char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= N - M; i++) {
        cilk_spawn NavieCal(pat,txt,i,M,N);
    }

    cilk_sync;
}



void PPNaive(char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i <= N - M; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < M; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == M) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            cout << "Pattern found at index " << i << endl;
    }
}


void Naive(char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= N - M; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < M; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == M) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            cout << "Pattern found at index " << i << endl;
    }
}


void computeLPSArray(char* pat, int M, int* lps);

// Prints occurrences of txt[] in pat[]
void KMPSearch(char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);

    // create lps[] that will hold the longest prefix suffix
    // values for pattern
    int* lps = new int[M];

    // Preprocess the pattern (calculate lps[] array)
    computeLPSArray(pat, M, lps);

    int i = 0; // index for txt[]
    int j = 0; // index for pat[]
    while ((N - i) >= (M - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }

        if (j == M) {
            cout << "Pattern found at index " << i - j << endl;
            j = lps[j - 1];
        }

        // mismatch after j matches
        else if (i < N && pat[j] != txt[i]) {
            // Do not match lps[0..lps[j-1]] characters,
            // they will match anyway
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
}

// Fills lps[] for given pattern pat[0..M-1]
void computeLPSArray(char* pat, int M, int* lps)
{
    // length of the previous longest prefix suffix
    int len = 0;

    lps[0] = 0; // lps[0] is always 0

    // the loop calculates lps[i] for i = 1 to M-1
    int i = 1;
    while (i < M) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        }
        else // (pat[i] != pat[len])
        {
            // This is tricky. Consider the example.
            // AAACAAAA and i = 7. The idea is similar
            // to search step.
            if (len != 0) {
                len = lps[len - 1];

                // Also, note that we do not increment
                // i here
            }
            else // if (len == 0)
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}




void RabinKarp(char* pat, char* txt, int q)
{
    int M = strlen(pat);
    int N = strlen(txt);
    int i, j;
    int p = 0;  // hash value for pattern
    int t = 0; // hash value for txt
    int h = 1;

    // The value of h would be "pow(d, M-1)%q"
    for (i = 0; i < M - 1; i++)
        h = (h * d) % q;

    // Calculate the hash value of pattern and first window of text
    for (i = 0; i < M; i++)
    {
        p = (d * p + pat[i]) % q;
        t = (d * t + txt[i]) % q;
    }

    // Slide the pattern over text one by one 
    for (i = 0; i <= N - M; i++)
    {

        // Chaeck the hash values of current window of text and pattern
        // If the hash values match then only check for characters on by one
        if (p == t)
        {
            /* Check for characters one by one */
            for (j = 0; j < M; j++)
            {
                if (txt[i + j] != pat[j])
                    break;
            }
            if (j == M)  // if p == t and pat[0...M-1] = txt[i, i+1, ...i+M-1]
            {
                printf("Pattern found at index %d \n", i);
            }
        }

        // Calulate hash value for next window of text: Remove leading digit, 
        // add trailing digit           
        if (i < N - M)
        {
            t = (d * (t - txt[i] * h) + txt[i + M]) % q;

            // We might get negative value of t, converting it to positive
            if (t < 0)
                t = (t + q);
        }
    }
}

#define REHASH(a, b, h) ((((h) - (a)*da) << 1) + (b))

int RabinKarp2(char* x, char* y) {
    int m = strlen(x);
    int n = strlen(y);
    int da, hx, hy, i, j, count;

    
        count = 0;
    /* Preprocessing */
    for (da = i = 1; i < m; ++i)
        da = (da << 1);

    for (hy = hx = i = 0; i < m; ++i) {
        hx = ((hx << 1) + x[i]);
        hy = ((hy << 1) + y[i]);
    }
    

    j = 0;
    while (j <= n - m) {
        if (hx == hy && memcmp(x, y + j, m) == 0) 
            cout << "Pattern found at index " << j << endl;
        hy = REHASH(y[j], y[j + m], hy);
        ++j;
    }
    
        return count;
}
void openfs() {

    std::cout << "Hello World!\n";
    string line;
    streampos size;
    char* memblock;
    ifstream myfile("bible.txt", ios::in | ios::binary | ios::ate);
    if (myfile.is_open())
    {
        size = myfile.tellg();
        memblock = new char[size];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, size);
        myfile.close();

        cout << "the entire file content is in memory";
        cout << memblock;

        delete[] memblock;
        /*
        while (getline(myfile, line))
        {
            cout << line << '\n';
        }
        myfile.close();*/
    }
    else cout << "Unable to open file";
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

int sol(char* x, char* y) {
    int m = strlen(x);
    int n = strlen(y);
    unsigned int lim, D, k, h, p_len;
    unsigned int S[256];
    int j, count;

    p_len = m;
    m = 32;

    /* Preprocessing */
    lim = preSo(x, m, S);

        /* Searching */
    count = 0;
    for (D = ~0, j = 0; j < n; ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            k = 0;
            h = j - m + 1;
            while (k < p_len && x[k] == y[h + k]) k++;
            if (k == p_len) 
                cout << "Pattern found at index " << j - m + 1 << endl;
        }
    }
    return count;
}





void preBmBc(char* x, int m, int bmBc[]) {
    int i;
    for (i = 0; i < SIGMA; ++i)
        bmBc[i] = m;
    for (i = 0; i < m - 1; ++i)
        bmBc[x[i]] = m - i - 1;
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


int bm(char* x, char* y) {
    int m = strlen(x);
    int n = strlen(y);
    int i, j, bmGs[XSIZE], bmBc[SIGMA], count;

    /* Preprocessing */
    preBmGs(x, m, bmGs);
    preBmBc(x, m, bmBc);
    

        /* Searching */
    j = 0;
    count = 0;
    while (j <= n - m) {
        for (i = m - 1; i >= 0 && x[i] == y[i + j]; --i);
        if (i < 0) {
            cout << "Pattern found at index " << j << endl;
            j += bmGs[0];
        }
        else
            j += max(bmGs[i], bmBc[y[i + j]] - m + 1 + i);
    }
    return count;
}



int so(char* x, char* y) {
    int m = strlen(x);
    int n = strlen(y);
    unsigned int lim, D;
    unsigned int S[256];
    int j, count;
    if (m > 32) return sol(x, y);

    /* Preprocessing */
    lim = preSo(x, m, S);

    count = 0;
    for (D = ~0, j = 0; j < n; ++j) {
        D = (D << 1) | S[y[j]];
        //cout << "D is " << D << endl;
        if (D < lim)
            cout << "Pattern found at index " << j - m + 1 << endl;
    }
    return count;
}
int soa(char* x, char* y) {
    int m = strlen(x);
    int n = strlen(y);
    unsigned int lim, D;
    unsigned int S[256];
    int j, count;
    if (m > 32) return sol(x, y);

    /* Preprocessing */
    lim = preSo(x, m, S);

    count = 0;
    D = ~0;

    for (j = 0; j < n; ++j) {
        D = (D << 1) | S[y[j]];
        //cout << "D is " << D << endl;
        if (D < lim)
            cout << "Pattern found at index " << j - m + 1 << endl;
    }
    return count;
}



void SNaive(char* pat, char* txt,int end)
{
    int M = strlen(pat);
    int N = strlen(txt);

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= end; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < M; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == M) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            cout << "Pattern found at index " << i << endl;
    }
}


void slice(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    for (int i = 0; i < tread; i++) {
        
        struct timeval t1, t2;
        gettimeofday(&t1,0);
        char* part = new char[(1) * (N / tread) + 1];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) );
        part[(1) * (N / tread)] = '\0';
        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        //printf("%f\n", runtime_ms/1000.0);
        //cout << i << endl;
        //cout << i * (N / tread) << endl;
        //cout << (i + 1) * (N / tread) << endl;
        cilk_spawn Naive(pat, part);
    }

    cilk_sync;
}


void slice2(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i < tread; i++) {
        
        struct timeval t1, t2;
        gettimeofday(&t1,0);
        char* part = new char[(1) * (N / tread) + 1];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) );
        part[(1) * (N / tread)] = '\0';
        //cout << i << endl;
        //cout << i * (N / tread) << endl;
        //cout << (i + 1) * (N / tread) << endl;
        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        //printf("%f\n", runtime_ms/1000.0);
        Naive(pat, part);
    }

    cilk_sync;
}




void sliceRabinKarp2(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i < tread; i++) {
        
        struct timeval t1, t2;
        gettimeofday(&t1,0);
        char* part = new char[(1) * (N / tread) + 1];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) );
        part[(1) * (N / tread)] = '\0';
        //cout << i << endl;
        //cout << i * (N / tread) << endl;
        //cout << (i + 1) * (N / tread) << endl;
        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        //printf("%f\n", runtime_ms/1000.0);
        RabinKarp2(pat, part);
    }

    cilk_sync;
}





void ShiftOrMatch(char* pSrc, char* pSubSrc)
{
    int nSrcSize = strlen(pSrc);
    int nSubSrcSize = strlen(pSubSrc);
    long skip[256];
    memset(skip, -1, sizeof(skip));
    for (int i = 0; i < nSubSrcSize; i++)
    {
        skip[pSubSrc[i]] ^= (0x01 << i);
    }

    long mask = ~(0x01 << (nSubSrcSize - 1));
    long ds = -1;
    int nPos = 0;
    while (nPos <= nSrcSize - nSubSrcSize)
    {
        ds = (ds << 0x01) | skip[pSrc[nPos]];
        if (~(ds | mask))
        {
            //break;
            cout << "Pattern found at index " << nPos << endl;
        }
        nPos++;
    }
    //return nPos - (nSubSrcSize - 1);

}

void sliceShiftOrMatch(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i < tread; i++) {
        
        struct timeval t1, t2;
        gettimeofday(&t1,0);
        char* part = new char[(1) * (N / tread) + 1];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) );
        part[(1) * (N / tread)] = '\0';
        //cout << i << endl;
        //cout << i * (N / tread) << endl;
        //cout << (i + 1) * (N / tread) << endl;
        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        //printf("%f\n", runtime_ms/1000.0);
        ShiftOrMatch(pat, part);
    }

    cilk_sync;
}
int main()
{
    
    std::cout << "Hello World!\n";
    string line;
    streampos size;
    char* memblock;
    ifstream myfile("text2.txt", ios::in | ios::binary | ios::ate);
    if (myfile.is_open())
    {
        size = myfile.tellg();
        memblock = new char[size];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, size);
        myfile.close();
        //char* text;
        //text = new char[size*10];
        //text = memblock + memblock + memblock + memblock + memblock + memblock + memblock + memblock + memblock + memblock;

        cout << "the entire file content is in memory"<< endl;
        //cout << memblock;

        struct timeval t1, t2;
        char pat[] = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababab";
        gettimeofday(&t1,0);
        //Naive(pat, memblock);
        gettimeofday(&t2,0);
        unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "Naive " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        PPNaive(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "PPNaive " << endl;
        printf("%f\n", runtime_ms/1000.0);

        gettimeofday(&t1,0);
        slice(pat, memblock,16);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "Slice Naive " << endl;
        printf("%f\n", runtime_ms/1000.0);


        gettimeofday(&t1,0);
        sliceRabinKarp2(pat, memblock,16);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << " sliceRabinKarp2 " << endl;
        printf("%f\n", runtime_ms/1000.0);


        gettimeofday(&t1,0);
        sliceShiftOrMatch(pat, memblock,16);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "sliceShiftOrMatch " << endl;
        printf("%f\n", runtime_ms/1000.0);


        gettimeofday(&t1,0);
        slice2(pat, memblock,16);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "Slice2 Naive " << endl;
        printf("%f\n", runtime_ms/1000.0);


        
        gettimeofday(&t1,0);
        //PNaive(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "PNaive " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        KMPSearch(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "KMPSearch " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        RabinKarp(pat, memblock, 7);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "RabinKarp " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        RabinKarp2(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "RabinKarp2 " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        //so(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "Shift Or" << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        //soa(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "Shift Or attemp" << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        bm(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "bm " << endl;
        printf("%f\n", runtime_ms/1000.0);
        gettimeofday(&t1,0);
        ShiftOrMatch(pat, memblock);
        gettimeofday(&t2,0);
        runtime_ms = (todval(&t2)-todval(&t1))/1000;
        cout << "ShiftOrMatch " << endl;
        printf("%f\n", runtime_ms/1000.0);


        delete[] memblock;
    }
    else cout << "Unable to open file";
    return 0;
}
