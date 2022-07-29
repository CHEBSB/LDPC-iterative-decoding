#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define iter_Max 100        // iteration limit of decoding
#define bSNR_dB 1.584          // Eb/N0 in dB

typedef struct edge {
    int v;      // the var node it connects
    int c;      // the check node it connects
    double q;   // bottom-up LLR
    double u;   // top-down LLR
} edge;

// seed for generating random number in (0, 1)
const unsigned long long SEED = 1024; 
unsigned long long RANV;
int RANI = 0;
double n1, n2;              // gaussian noise
// standard deviation of Gaussian noise
const double std = pow(10, bSNR_dB / ((double)-20)); 
// for LDPC
int n, r;                   // # of column/row in H     
int dv, dc;                 // # of 1 in column/row
int **col;                  // 1's in each column
int **row;                  // 1's in each row
double *L;
double *q;
int *x_hat;                 // decoder's output
edge*** V;          // edges from viewpoint of var node
edge*** C;          // edges from viewpoint of check node
int **A;                    /* for Bottom-up
A[i] = {0, 1, ..., i - 1, i + 1, i + 2, ..., dc - 1} */

// generate a uniform random number
double Ranq1();            
// generate 2 normal random number given standard deviation
void normal();
// value of H[i][j] (H is the parity-check matrix)
int H(int i, int j);
// find the edge connecting ci and vj
edge* Edge(int i, int j);
// check node operation
double CHK(double L1, double L2);
// iterative decoding algorithm
int iterDecod(double *y);

int main(void)
{
    int i, j, k;            // looping indices
    int temp;               // temp storage
    // previous information bits
    // int U[] = {0, 0, 0, 0, 0, 0};
    // int u;                  // current information bit
    // int PN[63];             // 1 period of the PN sequence
    // double *x;           // encoder output
    double *y;              // decoder's input sequence
    int err = 0;            // # of errors
    int err_block = 0;      // # of blocks with error
    int decode;             // boolean for decoder failure

    // read input from file
    scanf("%d %d", &n, &r);
    scanf("%d %d", &dv, &dc);
    // allocate memory
    A = (int **)calloc(dc, sizeof(int *));
    for (i = 0; i < dc; i++) 
        A[i] = (int *)calloc(dc - 1, sizeof(int));
    y = (double *)calloc(n, sizeof(double));
    col = (int **)calloc(n, sizeof(int *));
    row = (int **)calloc(r, sizeof(int *));
    V = (edge ***)calloc(n, sizeof(edge **));
    C = (edge ***)calloc(r, sizeof(edge **));
    L = (double *)calloc(n, sizeof(double));
    q = (double *)calloc(n, sizeof(double));
    x_hat = (int *)calloc(n, sizeof(int));
    for (i = 0; i < r; i++) {
        row[i] = (int *)calloc(dc, sizeof(int));
        C[i] = (edge **)calloc(dc, sizeof(edge *));
    }
    for (i = 0; i < n; i++) {
        col[i] = (int *)calloc(dv, sizeof(int));
        V[i] = (edge **)calloc(dv, sizeof(edge *));
    }
    // skip the following (n + r) numbers
    for (i = 0; i < n + r; i++)
        scanf("%d", &j);
    // for each column
    for (i = 0; i < n; i++)
        for (j = 0; j < dv; j++) {
            scanf("%d", &temp);
            col[i][j] = temp - 1;
            V[i][j] = (edge *)malloc(sizeof(edge));
            V[i][j]->v = i;
            V[i][j]->c = temp - 1;
        }
    // for each row
    for (i = 0; i < r; i++) 
        for (j = 0; j < dc; j++) {
            scanf("%d", &temp);
            row[i][j] = temp - 1;
            C[i][j] = Edge(i, temp - 1);
        }
    // init A
    for (i = 0; i < dc; i++) {
        k = 0;
        for (j = 0; j < dc; j++)
            if (i != j) {
                A[i][k] = j;
                k++;
            }
    }
    printf("n = %d r = %d dv = %d dc = %d\n", n, r, dv, dc);   // for debug
    // test until # of blocks with error reach 50
    for (i = 0; err_block < 50; i++) {
        // generate decoder's input
        for (j = 0; j < n; j += 2) {
            normal();           // generate Gaussian noise
            y[j] = 1 + n1;      // all-zero input
            if (j + 1 < n) y[j + 1] = 1 + n2;   
        }
        decode = iterDecod(y);
        if (decode == 1) {
            // err_block += 1;
            // err += (n - r);
            printf("%d-th Failure at %d-th!\n", err_block + 1, i);
        }
        //else {
        temp = 0;
        for (j = 0; j < n - r; j++)
            if (x_hat[j] != 0) {
                err++;
                temp = 1;
            }
        err_block += temp;
        //}
        if (i % 100 == 0) printf("i = %d\terr: %d\n", i, err);
    }
    printf("BER = %lf * 10^-3\n", ((double)err) / (n - r) * 1000 / (i + 1));
    return 0;
} 

// generate an uniform random number in (0, 1)
double Ranq1()
{
    if (RANI == 0) {
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

// generate 2 normal random number given standard deviation
void normal()
{
    double x1, x2, s;       // x1 and x2 are Uniform(0, 1)

    do {
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;    // transform (0, 1) to (-1, 1)
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while (s >= 1.0);     // until s is in unit circle
    n1 = std * x1 * sqrt(-2 * log(s) / s);
    n2 = std * x2 * sqrt(-2 * log(s) / s);
}

// value of H[i][j] (H is the parity-check matrix)
int H(int i, int j)
{
    int k;              // looping index

    // check if there exist k s.t. col[j][k] = i
    for (k = 0; col[j][k] != i && k < dv; k++);
    if (k == dv) return 0;      // not found
    return 1;
}

// find the edge connecting ci and vj
edge* Edge(int i, int j)
{
    int k;          // looping index

    for (k = 0; k < dv; k++) 
        if (V[j][k]->c == i) 
            return V[j][k];
    return NULL;
}

// check node operation
double CHK(double L1, double L2)
{
    int s1, s2;     // sign of L1, L2
    double A1, A2;  // A1 = |L1|, A2 = |L2|
    double delta;   // the complicated term

    delta = log((1 + exp(-1 * fabs(L1 + L2))) 
    / (1 + exp(-1 * fabs(L1 - L2))));
    A1 = fabs(L1);
    A2 = fabs(L2);
    s1 = (L1 >= 0)? 1: -1;
    s2 = (L2 >= 0)? 1: -1;
    if (A1 > A2)
        return s1 * s2 * A2 + delta;
    return s1 * s2 * A1 + delta;
}

// iterative decoding
int iterDecod(double *y)
{

    int i, j, k, m;         // looping indices
    int check;              // value for parity check
    int flag;               // if parity-check continues

    // Lj = 2 * yj / (std ^ 2)
    for (j = 0; j < n; j++) 
        L[j] = 2 * y[j] / std / std;
    // init: for each edge, set qij = Lj
    for (j = 0; j < n; j++) 
        for (i = 0; i < dv; i++)
            V[j][i]->q = L[j];
    // at most 100 iteration
    for (k = 0; k < iter_Max; k++) {
    // Bottom-up
        // for each check node ci
        for (i = 0; i < r; i++) 
            // for each edge connected with ci
            for (j = 0; j < dc; j++) {
                C[i][j]->u = CHK(C[i][A[j][0]]->q, C[i][A[j][1]]->q);
                for (m = 2; m < dc - 1; m++)
                    C[i][j]->u = CHK(C[i][j]->u, C[i][A[j][m]]->q);
            }
    // Top-down
        // for each var node vj
        for (j = 0; j < n; j++)
            // for each edge connected with vj
            for (i = 0; i < dv; i++) {
                V[j][i]->q = L[j];
                for (m = 0; m < dv; m++) 
                    if (m != i) V[j][i]->q += V[j][m]->u;
            }
    // Decision
        // compute qj
        for (j = 0; j < n; j++) {
            q[j] = L[j];
            for (i = 0; i < dv; i++) {
                q[j] += V[j][i]->u;
            }
            if (q[j] >= 0)
                x_hat[j] = 0;
            else
                x_hat[j] = 1;
        }
        // check (x_hat) by each parity-check equation
        flag = 1;
        for (i = 0; i < r && flag == 1; i++) {
            check = 0;
            for (j = 0; j < dc; j++) 
                check += x_hat[row[i][j]];
            if ((check % 2) != 0) 
                flag = 0;
        }
        if (flag == 1)          // x_hat be the output
            return 0;           // terminate
    }
    return 1;       // decoding failure
}