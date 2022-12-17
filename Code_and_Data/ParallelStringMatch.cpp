#include <iostream>
#include <string.h>
#include <fstream>

#define d 256 
#define XSIZE 4200
#define SIGMA 256

using namespace std;

#include <sys/time.h>
#include <cilk/cilk.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <random>

#define REHASH(a, b, h) ((((h) - (a)*da) << 1) + (b))


unsigned long long todval (struct timeval *tp) {
    return tp->tv_sec * 1000 * 1000 + tp->tv_usec;
}




int PNaive(char* pat, char* txt,int patlen, int txtlen)
{
    int co = 0;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i <= txtlen - patlen; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < patlen; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == patlen) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            {
                cout << i;
                co++;
            }

    }
    return co;
}


int Naive(char* pat, char* txt,int patlen, int txtlen)
{
    int co = 0;
    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= txtlen - patlen; i++) {
        int j;

        /* For current index i, check for pattern match */
        for (j = 0; j < patlen; j++)
            if (txt[i + j] != pat[j])
                break;

        if (j == patlen) // if pat[0...M-1] = txt[i, i+1, ...i+M-1]
            {
                cout << i;
                co++;
            }
    }
    return co;
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

int KMPSearch(char* pat, char* txt,int patlen, int txtlen)
{
    int* lps = new int[patlen];
    int co = 0;
    LPSArray(pat, patlen, lps);

    int i = 0; 
    int j = 0; 
    while ((txtlen - i) >= (patlen - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == patlen) {
            cout <<  i - j;
            j = lps[j - 1];
            co++;
        }
        else if (i < txtlen && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
    return co;
}

int KMPSearch(char* pat, char* txt, int patlen, int txtlen, int* lps)
{
    int co = 0;
    int i = 0;
    int j = 0;
    while ((txtlen - i) >= (patlen - j)) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == patlen) {
            cout << i - j;
            j = lps[j - 1];
            co++;
        }
        else if (i < txtlen && pat[j] != txt[i]) {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
    return co;
}

int RabinKarp2(char* x, char* y,int patlen, int txtlen)
{
    int da, hx, hy, i, j, count;
    count = 0;
    for (da = i = 1; i < patlen; ++i)
        da = (da << 1);

    for (hy = hx = i = 0; i < patlen; ++i) {
        hx = ((hx << 1) + x[i]);
        hy = ((hy << 1) + y[i]);
    }

    

    j = 0;
    while (j <= txtlen - patlen) {
        if (hx == hy && memcmp(x, y + j, patlen) == 0) {

            cout << j;
            count++;
        }
        hy = REHASH(y[j], y[j + patlen], hy);
        ++j;
    }
    
        return count;
}


void openfs(const char* name,char* memblock) {
    string line;
    streampos size;
    ifstream myfile(name, ios::in | ios::binary | ios::ate);
    if (myfile.is_open())
    {
        size = myfile.tellg();
        memblock = new char[size];
        myfile.seekg(0, ios::beg);
        myfile.read(memblock, size);
        myfile.close();

        cout << "the entire file content is in memory";
        //cout << memblock;

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


int sol(char* x, char* y,int patlen, int txtlen)
{
    unsigned int lim, D, k, h, p_len;
    unsigned int S[256];
    int j, count;

    p_len = patlen;
    patlen = 32;

    /* Preprocessing */
    lim = preSo(x, patlen, S);

        /* Searching */
    count = 0;
    for (D = ~0, j = 0; j < txtlen; ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim) {
            k = 0;
            h = j - patlen + 1;
            while (k < p_len && x[k] == y[h + k]) k++;
            if (k == p_len) {
                count++;
                cout <<  j - patlen + 1 ;
            }
        }
    }
    return count;
}

int so(char* x, char* y,int patlen, int txtlen)
{
    unsigned int lim, D;
    unsigned int S[256];
    int j, count;
    if (patlen > 32) return sol(x, y, patlen,  txtlen);
    lim = preSo(x, patlen, S);
    count = 0;
    for (D = ~0, j = 0; j < txtlen; ++j) {
        D = (D << 1) | S[y[j]];
        if (D < lim){
            count++;
            cout << j - patlen + 1 ;
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


int bm(char* x, char* y,int patlen, int txtlen)
{
    int i, j, bmGs[XSIZE], bmBc[SIGMA], count;

    /* Preprocessing */
    preBmGs(x, patlen, bmGs);
    preBmBc(x, patlen, bmBc);
    

        /* Searching */
    j = 0;
    count = 0;
    while (j <= txtlen - patlen) {
        for (i = patlen - 1; i >= 0 && x[i] == y[i + j]; --i);
        if (i < 0) {
            cout <<  j ;
            j += bmGs[0];
        }
        else
            j += max(bmGs[i], bmBc[y[i + j]] - patlen + 1 + i);
    }
    return count;
}



void sliceNaive(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M ];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        Naive(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }
    cilk_sync;
}

void sliceKMP(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    int* lps = new int[strlen(pat)];
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M ];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        KMPSearch(pat, part, M, (1) * (N / tread) + M , lps);
        delete[] part;
    }

    cilk_sync;
}


void sliceRabinKarp2(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M ];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        RabinKarp2(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }

    cilk_sync;
}



void sliceBM(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M ];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        bm(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }
    cilk_sync;
}




void sliceso(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M ];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        so(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }
    cilk_sync;
}


typedef union {
    __m128i* data16;
    unsigned char* data;
} TEXT;

typedef struct list {
    struct list* next;
    int pos;
} LIST;


int SSEF(char* x, char* y,int Plen,int Tlen) {
    if (Plen < 32) return -1;
    LIST* flist[65536];
    LIST* t;
    memset(flist, 0, sizeof(LIST*) * 65536);
    __m128i tmp128;
    TEXT T;
    T.data16 = (__m128i*) y;
    T.data = (unsigned char*)y;

    
    unsigned int K = 7;
    unsigned int count = 0;
    unsigned int i, last, j;
    __m128i* ptr16;
    __m128i* lastchunk = &T.data16[Tlen / 16];
    unsigned int filter;
    unsigned char* f = (unsigned char*)malloc(Plen);

    last = (Plen / 16) - 1;
    for (i = 0; i < Plen; i++) {
        f[i] = (x[i] & 0x80) >> 7;
    }
    count = 15;

    for (i = 0; i < last * 16; i++) {
        j = last * 16 - i;
        filter = f[j] + f[j + 1] * 2 + f[j + 2] * 4 + f[j + 3] * 8 +
            f[j + 4] * 16 + f[j + 5] * 32 + f[j + 6] * 64 + f[j + 7] * 128 +
            f[j + 8] * 256 + f[j + 9] * 512 + f[j + 10] * 1024 + f[j + 11] * 2048 +
            f[j + 12] * 4096 + f[j + 13] * 8192 + f[j + 14] * 16384 + f[j + 15] * 32768;
        if (flist[filter] == 0) {
            flist[filter] = (LIST*)malloc(sizeof(LIST));
            flist[filter]->next = NULL;
            flist[filter]->pos = i;
        }
        else {
            t = flist[filter];
            while (t->next != NULL) t = t->next;
            t->next = (LIST*)malloc(sizeof(LIST));
            t = t->next;
            t->next = NULL;
            t->pos = i;
        }
    }
    

        
    count = 0;
    ptr16 = &T.data16[last];

    while (ptr16 < lastchunk) {
        filter = _mm_movemask_epi8(*ptr16);

        if (flist[filter]) {
            i = ((ptr16 - &T.data16[0]) - last) * 16;
            t = flist[filter];
            while (t) {
                if (memcmp(x, &T.data[i + t->pos], Plen) == 0) {
                    count++;
                    cout <<last;

                }
                t = t->next;
            }
        }
        ptr16 += last;
    }
    return count;
}




typedef union {
    __m128i  v;
    unsigned  int  ui[4];
    unsigned short int  us[8];
    unsigned char  uc[16];
}VectorUnion;

int EPSM1(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{//we exactly know patlen=1

    __m128i* text = (__m128i*)x;
    __m128i* end = (__m128i*)(x + 16 * (textlen / 16));
    __m128i t0, a;
    VectorUnion template0;
    unsigned int j, k;
    int cnt = 0;


    for (j = 0; j < 16; j++)
    {
        template0.uc[j] = pattern[0];
    }
    t0 = template0.v;



    while (text < end) {
        a = _mm_cmpeq_epi8(t0, *text);
        j = _mm_movemask_epi8(a);
        cnt += _mm_popcnt_u32(j);
        text++;
    }
    //now we are at the beginning of the last 16-byte block, perform naive check
    for (j = 16 * (textlen / 16); j < textlen; j++)
        cnt += (x[j] == pattern[0]);

    return cnt;
}

int EPSM2(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{//we exactly know patlen=2

    __m128i* text = (__m128i*)x;
    __m128i* end = (__m128i*)(x + 16 * (textlen / 16));
    __m128i t0, t1, a, b;
    VectorUnion template0, template1;
    unsigned int j, k, carry = 0;
    int cnt = 0;
    unsigned char firstch = pattern[0], lastch = pattern[1];


    for (j = 0; j < 16; j++)
    {
        template0.uc[j] = firstch;        //template0.uc[i+1]=lastch;
        template1.uc[j] = lastch;        //template1.uc[i+1]=firstch;
    }
    t0 = template0.v;
    t1 = template1.v;



    while (text < end) {
        a = _mm_cmpeq_epi8(t0, *text);
        j = _mm_movemask_epi8(a);
        b = _mm_cmpeq_epi8(t1, *text);
        k = _mm_movemask_epi8(b);
        cnt += _mm_popcnt_u32(((j << 1) | (carry >> 15)) & k);
        carry = j & 0x00008000;
        text++;
    }
    //now we are at the beginning of the last 16-byte block, perform naive check
    for (j = 16 * (textlen / 16); j < textlen; j++)
        cnt += ((x[j - 1] == firstch) && (x[j] == lastch));

    return cnt;
}

int EPSM3(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{//we exactly know patlen=3


    __m128i* text = (__m128i*)x;
    __m128i* end = (__m128i*)(x + 16 * (textlen / 16));
    __m128i t0, t1, t2, a, b, c;
    VectorUnion template0, template1, template2;
    unsigned int j, k, l, carry0 = 0, carry1 = 0;
    int cnt = 0;
    for (j = 0; j < 16; j++)
    {
        template0.uc[j] = pattern[0];
        template1.uc[j] = pattern[1];
        template2.uc[j] = pattern[2];
    }
    t0 = template0.v;
    t1 = template1.v;
    t2 = template2.v;



    while (text < end) {
        a = _mm_cmpeq_epi8(t0, *text);
        j = _mm_movemask_epi8(a);

        b = _mm_cmpeq_epi8(t1, *text);
        k = _mm_movemask_epi8(b);

        c = _mm_cmpeq_epi8(t2, *text);
        l = _mm_movemask_epi8(c);

        cnt += _mm_popcnt_u32(((j << 2) | (carry0 >> 14)) & ((k << 1) | (carry1 >> 15)) & l);
        carry0 = j & 0x0000C000;
        carry1 = k & 0x00008000;
        text++;
    }
    //now we are at the beginning of the last 16-byte block, perform naive check
    for (j = 16 * (textlen / 16); j < textlen; j++)
        cnt += ((x[j - 2] == pattern[0]) && (x[j - 1] == pattern[1]) && (x[j] == pattern[2]));

    return cnt;
}

int EPSM4(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    __m128i* text = (__m128i*)x;
    __m128i* end = (__m128i*)(x + 16 * (textlen / 16));
    if ((textlen % 16) < 7) end--;


    int i, count = 0;
    VectorUnion P, Z;
    __m128i a, b, p, z;

    Z.ui[0] = Z.ui[1] = Z.ui[2] = Z.ui[3] = 0;
    z = Z.v;
    P.uc[0] = pattern[0];
    P.uc[1] = pattern[1];
    P.uc[2] = pattern[2];
    P.uc[3] = pattern[3];
    p = P.v;

    text++;// leave the naive check of the first block to the end



    while (text < end)
    {
        //check if P[(m-4) ... (m-1)] matches with T[i*16 ... i*16+3], T[i*16+1 ... i*16+4], .... , T[i*16+7 ... i*16+10]
        a = _mm_mpsadbw_epu8(*text, p, 0x00);
        b = _mm_cmpeq_epi16(a, z);
        i = _mm_movemask_epi8(b);
        count += _mm_popcnt_u32(i);

        a = _mm_blend_epi16(*text, *(text + 1), 0x0f);
        b = _mm_shuffle_epi32(a, _MM_SHUFFLE(1, 0, 3, 2));

        //check if P[(m-4) ... (m-1)] matches with T[i*16+8 ... i*16+11], T[i*16+9 ... i*16+12], .... , T[i*16+15 ... i*16+18]
        a = _mm_mpsadbw_epu8(b, p, 0x00);
        b = _mm_cmpeq_epi16(a, z);
        i = _mm_movemask_epi8(b);
        count += _mm_popcnt_u32(i);
        text++;
    }
    count = count / 2;

    // the ending position of the pattern from the first appropriate position T[patlen-1] to the third position of the next 16-byte block is performed naive
    for (i = 3; (i < 19) && (i < textlen); i++) // j presents possible end points of the pattern
        if (0 == memcmp(pattern, &x[i - 3], patlen)) count++;

    // note that at the last iteration of the while loop, we have checked if P ends at positions 0,1,and 2 of the last 16-byte block
    // however, what if the last position of the text is beyond 2?
    // for the possibilities that T ends at positions 3,4,5,6,7,8,9,10,11,12,13,14,and 15, we perform naive checks

    for (i = ((unsigned char*)text) + 3 - x; i < textlen; i++)
    {
        if (0 == memcmp(pattern, &x[i - 3], 4)) count++;
    }


    return count;
}

int EPSM16(unsigned char* pattern, int patlen, unsigned char* x, int textlen)
{
    LIST* flist[2048]; //11 bit hash is gives the best result according to our tests, no shorter no longer
    LIST* t, * t1;


    unsigned int i, filter, shift = (patlen / 8) - 1;
    unsigned long long crc, seed = 123456789, mask;
    unsigned long long* ptr64;
    unsigned long long* lastchunk;
    unsigned char* charPtr;
    int count = 0, p = 0;
    mask = 2047;

    memset(flist, 0, sizeof(LIST*) * 2048);

    for (i = 1; i < patlen - 7; i++)
    {
        ptr64 = (unsigned long long*)(&pattern[i]);
        crc = _mm_crc32_u64(seed, *ptr64);
        filter = (unsigned int)(crc & mask);

        if (flist[filter] == 0)
        {
            flist[filter] = (LIST*)malloc(sizeof(LIST));
            if (flist[filter]) {
                flist[filter]->next = 0;
                flist[filter]->pos = i;
            }
        }
        else
        {
            t = flist[filter];
            while (t->next != 0) t = t->next;
            t->next = (LIST*)malloc(sizeof(LIST));
            if (t->next) {
                t = t->next;
                t->next = 0;
                t->pos = i;
            }
        }
    }

    lastchunk = (unsigned long long*) & x[((textlen - patlen) / 8) * 8 + 1];
    ptr64 = (unsigned long long*) & x[(shift - 1) * 8];

    crc = _mm_crc32_u64(seed, *ptr64);
    filter = (unsigned int)(crc & mask);
    if (flist[filter]) {
        charPtr = (unsigned char*)ptr64;
        t = flist[filter];
        while (t)
        {

            if (t->pos <= 8 * (shift - 1)) {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
            }
            t = t->next;
        }
    }
    ptr64 += shift;

    crc = _mm_crc32_u64(seed, *ptr64);
    filter = (unsigned int)(crc & mask);
    if (flist[filter]) {
        charPtr = (unsigned char*)ptr64;
        t = flist[filter];
        while (t)
        {

            if (t->pos <= 8 * (2 * shift - 1)) {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
            }
            t = t->next;
        }
    }
    ptr64 += shift;



    while (ptr64 < lastchunk)
    {
        crc = _mm_crc32_u64(seed, *ptr64);
        filter = (unsigned int)(crc & mask);

        if (flist[filter])
        {
            charPtr = (unsigned char*)ptr64;
            t = flist[filter];
            while (t)
            {
                if (memcmp(pattern, charPtr - t->pos, patlen) == 0)
                    count++;
                t = t->next;
            }
        }
        ptr64 += shift;
    }

    ptr64 -= shift;
    charPtr = (unsigned char*)ptr64;
    charPtr += patlen - 1; //the first position unchecked where P may end

    while (charPtr < &x[textlen - 1])
    {
        if (0 == memcmp(pattern, charPtr - patlen + 1, patlen)) count++;
        charPtr++;
    }

    return count;
}


int EPSM(char* pattern, char* x,int patlen, int textlen)
{
    if (patlen < 2)         return EPSM1((unsigned char*)pattern, patlen, (unsigned char*)x, textlen);
    if (patlen == 2)        return EPSM2((unsigned char*)pattern, patlen, (unsigned char*)x, textlen);
    if (patlen == 3)        return EPSM3((unsigned char*)pattern, patlen, (unsigned char*)x, textlen);
    if (patlen == 4)        return EPSM4((unsigned char*)pattern, patlen, (unsigned char*)x, textlen);
    if (patlen >= 16)       return EPSM16((unsigned char*)pattern, patlen, (unsigned char*)x, textlen);


    unsigned char* y0;
    int i, j, k, count = 0;
    VectorUnion P, zero;
    __m128i res, a, b, z, p;

    __m128i* text = (__m128i*)x;
    __m128i* end = (__m128i*)(x + 16 * (textlen / 16));
    end--;

    zero.ui[0] = zero.ui[1] = zero.ui[2] = zero.ui[3] = 0;  z = zero.v;
    P.uc[0] = pattern[patlen - 5];  P.uc[1] = pattern[patlen - 4];  P.uc[2] = pattern[patlen - 3];  P.uc[3] = pattern[patlen - 2];  p = P.v;

    i = (patlen - 1) / 16; // i points the first 16-byte block that P may end in
    i++;
    text += i;
    for (k = 0; k < (i * 16 + 8) - patlen + 1; k++) 	if (0 == memcmp(pattern, x + k, patlen)) count++;



    //the loop checks if pattern ends at the second half of text[i] or at the first half of text[i+1]
    while (text < end)
    {
        //check if P[(m-5) ... (m-2)] matches with T[i*16+4 ... i*16+7], T[i*16+5 ... i*16+8], .... , T[i*16+11 ... i*16+14]
        //note thet this corresponds P ends at T[i*16+8],T[i*16+9],...,T[i*16+15]
        res = _mm_mpsadbw_epu8(*text, p, 0x04);
        b = _mm_cmpeq_epi16(res, z);
        j = _mm_movemask_epi8(b);
        if (j)
        {
            y0 = (unsigned char*)(text)+9 - patlen;
            if ((j & 3) == 3 && !memcmp(pattern, y0, patlen)) count++;
            if ((j & 12) == 12 && !memcmp(pattern, y0 + 1, patlen)) count++;
            if ((j & 48) == 48 && !memcmp(pattern, y0 + 2, patlen)) count++;
            if ((j & 192) == 192 && !memcmp(pattern, y0 + 3, patlen)) count++;
            if ((j & 768) == 768 && !memcmp(pattern, y0 + 4, patlen)) count++;
            if ((j & 3072) == 3072 && !memcmp(pattern, y0 + 5, patlen)) count++;
            if ((j & 12288) == 12288 && !memcmp(pattern, y0 + 6, patlen)) count++;
            if ((j & 49152) == 49152 && !memcmp(pattern, y0 + 7, patlen)) count++;
        }

        a = _mm_blend_epi16(*text, *(text + 1), 0x0f);
        b = _mm_shuffle_epi32(a, _MM_SHUFFLE(1, 0, 3, 2));

        //check if P ends at T[(i+1)*16+8],T[(i+1)*16+9],...,T[(i+1)*16+15]
        res = _mm_mpsadbw_epu8(b, p, 0x04);
        b = _mm_cmpeq_epi16(res, z);
        j = _mm_movemask_epi8(b);

        if (j)
        {
            y0 = (unsigned char*)(text)+9 + 8 - patlen;
            if ((j & 3) == 3 && !memcmp(pattern, y0, patlen)) count++;
            if ((j & 12) == 12 && !memcmp(pattern, y0 + 1, patlen)) count++;
            if ((j & 48) == 48 && !memcmp(pattern, y0 + 2, patlen)) count++;
            if ((j & 192) == 192 && !memcmp(pattern, y0 + 3, patlen)) count++;
            if ((j & 768) == 768 && !memcmp(pattern, y0 + 4, patlen)) count++;
            if ((j & 3072) == 3072 && !memcmp(pattern, y0 + 5, patlen)) count++;
            if ((j & 12288) == 12288 && !memcmp(pattern, y0 + 6, patlen)) count++;
            if ((j & 49152) == 49152 && !memcmp(pattern, y0 + 7, patlen)) count++;
        }
        text++;
    }

    for(k = ((char*)text)+8-x ; k < textlen ; k++)
    {
        if (0 == memcmp(pattern, &x[k - patlen + 1], patlen)) count++;
    }


    return count;
}



void sliceEPSM(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        EPSM(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }

    cilk_sync;
}

void sliceSSEF(char* pat, char* txt, int tread) {

    int M = strlen(pat);
    int N = strlen(txt);
    //cout << N << endl;

    /* A loop to slide pat[] one by one */
    cilk_for (int i = 0; i < tread; i++) {
        char* part = new char[(1) * (N / tread) + M];
        memcpy(part, txt + (i * (N / tread)), (1) * (N / tread) + M - 1);
        part[(1) * (N / tread) + M - 1] = '\0';
        SSEF(pat, part, M, (1) * (N / tread) + M);
        delete[] part;
    }

    cilk_sync;
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

    char* pat = new char[patsize+1];
    memcpy(pat, memblock + 1000, patsize);
    pat[patsize] = '\0';
    /*if (myfile3.is_open()){

        size = myfile3.tellg();
        pat = new char[size];
        myfile3.seekg(0, ios::beg);
        myfile3.read(pat, size);
        myfile3.close();
        cout << pat << endl;
    }*/
    cout << "pat gen"<< endl;
    ofstream outfile;
    outfile.open(outf);
    outfile << pat << endl;

    struct timeval t1, t2;

    int patlen = strlen(pat);
    int txtlen = strlen(memblock);
    gettimeofday(&t1,0);
    Loop(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    unsigned long long runtime_ms = (todval(&t2)-todval(&t1))/1000;
    cout << "Loop " << endl;
    printf("%f\n", runtime_ms/1000.0);


    gettimeofday(&t1,0);
    Naive(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "Naive: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceNaive(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "Slice Naive: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    PNaive(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "Parallel Naive: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    bm(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "bm: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceBM(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "Slice BM: " << runtime_ms/1000.0 << endl ;


    gettimeofday(&t1,0);
    KMPSearch(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "KMPSearch: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceKMP(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "sliceKMP: " << runtime_ms/1000.0 << endl ;


    gettimeofday(&t1,0);
    so(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "Shift Or: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceso(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "slice Shift Or: " << runtime_ms/1000.0 << endl ;


    gettimeofday(&t1,0);
    RabinKarp2(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "RabinKarp: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceRabinKarp2(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "sliceRabinKarp2: " << runtime_ms/1000.0 << endl ;


    gettimeofday(&t1,0);
    EPSM(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "EPSM: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceEPSM(pat, memblock,16);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "sliceEPSM: " << runtime_ms/1000.0 << endl ;


    gettimeofday(&t1,0);
    SSEF(pat, memblock,patlen,txtlen);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "SSEF: " << runtime_ms/1000.0 << endl ;
    gettimeofday(&t1,0);
    sliceSSEF(pat, memblock,n);
    gettimeofday(&t2,0);
    runtime_ms = (todval(&t2)-todval(&t1))/1000;
    outfile << endl << "sliceSSEF: " << runtime_ms/1000.0 << endl ;

        
    delete[] memblock;
}

int main(){
    does(2,(char*)"gen50.txt",2, (char*)"ogen50p2t2.txt");
    does(4,(char*)"gen50.txt",2, (char*)"ogen50p2t4.txt");
    does(8,(char*)"gen50.txt",2, (char*)"ogen50p2t8.txt");
    does(16,(char*)"gen50.txt",2, (char*)"ogen50p2t16.txt");
    does(32,(char*)"gen50.txt",2, (char*)"ogen50p2t32.txt");
    does(2,(char*)"gen50.txt",4, (char*)"ogen50p4t2.txt");
    does(4,(char*)"gen50.txt",4, (char*)"ogen50p4t4.txt");
    does(8,(char*)"gen50.txt",4, (char*)"ogen50p4t8.txt");
    does(16,(char*)"gen50.txt",4, (char*)"ogen50p4t16.txt");
    does(32,(char*)"gen50.txt",4, (char*)"ogen50p4t32.txt");
    does(2,(char*)"gen50.txt",8, (char*)"ogen50p8t2.txt");
    does(4,(char*)"gen50.txt",8, (char*)"ogen50p8t4.txt");
    does(8,(char*)"gen50.txt",8, (char*)"ogen50p8t8.txt");
    does(16,(char*)"gen50.txt",8, (char*)"ogen50p8t16.txt");
    does(32,(char*)"gen50.txt",8, (char*)"ogen50p8t32.txt");
    does(2,(char*)"gen50.txt",16, (char*)"ogen50p16t2.txt");
    does(4,(char*)"gen50.txt",16, (char*)"ogen50p16t4.txt");
    does(8,(char*)"gen50.txt",16, (char*)"ogen50p16t8.txt");
    does(16,(char*)"gen50.txt",16, (char*)"ogen50p16t16.txt");
    does(32,(char*)"gen50.txt",16, (char*)"ogen50p16t32.txt");
    does(2,(char*)"gen50.txt",32, (char*)"ogen50p32t2.txt");
    does(4,(char*)"gen50.txt",32, (char*)"ogen50p32t4.txt");
    does(8,(char*)"gen50.txt",32, (char*)"ogen50p32t8.txt");
    does(16,(char*)"gen50.txt",32, (char*)"ogen50p32t16.txt");
    does(32,(char*)"gen50.txt",32, (char*)"ogen50p32t32.txt");
    does(2,(char*)"gen50.txt",64, (char*)"ogen50p64t2.txt");
    does(4,(char*)"gen50.txt",64, (char*)"ogen50p64t4.txt");
    does(8,(char*)"gen50.txt",64, (char*)"ogen50p64t8.txt");
    does(16,(char*)"gen50.txt",64, (char*)"ogen50p64t16.txt");
    does(32,(char*)"gen50.txt",64, (char*)"ogen50p64t32.txt");
}