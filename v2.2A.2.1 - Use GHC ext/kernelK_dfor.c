#include "mex.h"
#include "stdint.h"
#include <math.h>

// --- DICHIARAZIONI FUNZIONI ---
static inline void nucleoKL(double (* restrict nuKL)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
static inline void nucleoKT(double (* restrict nuKT)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
static inline void nucleoKLT(double (* restrict nuKLT)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double* restrict nodeInt, const double (* restrict matCoeff)[3], const double* restrict vetCoeff, const size_t indSM, const double cP, const double cS, const double lambda, const double mu);
static inline void nucleoKRj(double (* restrict nuKRj)[3], const double* restrict x, const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const size_t j);
static inline void nucleoKR(double (* restrict nuKR)[3], const double* restrict x, const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const double* restrict normInt, const double (* restrict matCoeff)[3], const double* restrict vettVMS);

static inline double dotProd3D(const double vettA[3], const double vettB[3]) {return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];}
static inline double norm2(const double* restrict v) {return sqrt(dotProd3D(v, v));}
static inline double baseFunctionSM(const double* restrict nodeInt, const double (* restrict matCoeff)[3], const double* restrict vetCoeff, const size_t indSM);

static inline double pow2(const double x) {return x * x;}
static inline double pow3(const double x) {return x * x * x;}
static inline double pow5(const double x) {const double x2 = pow2(x); return x2 * x2 * x;}

// Permutation tensor
static const double tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, 
                                        {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, 
                                        {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};
// --- MEX FUNCTION ---

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs == 0)
        return;

    double* temp = 0x0;
    
    const size_t numNodes = (size_t) mxGetScalar(prhs[0]);

    const double* const intNodes = mxGetPr(prhs[1]);
    const double* const intWeights = mxGetPr(prhs[2]);
    
    temp = mxGetPr(prhs[3]);
    const double nodeExt[3] = {temp[0], temp[1], temp[2]};

    const double currT = mxGetScalar(prhs[4]);

    const double cP = mxGetScalar(prhs[5]);
    const double cS = mxGetScalar(prhs[6]);

    const double lambda = mxGetScalar(prhs[7]);
    const double mu = mxGetScalar(prhs[8]);
    const double rho = mxGetScalar(prhs[9]);
    
    temp = mxGetPr(prhs[10]);
    const double vettVMS[3] = {temp[0], temp[1], temp[2]};

    temp = mxGetPr(prhs[11]);
    const double matCoeff[3][3] = {{temp[0], temp[3], temp[6]}, 
                                   {temp[1], temp[4], temp[7]},
                                   {temp[2], temp[5], temp[8]}};

    temp = mxGetPr(prhs[12]);
    const double vetCoeff[3] = {temp[0], temp[1], temp[2]};

    const size_t indSM = (size_t) mxGetScalar(prhs[13]);

    temp = mxGetPr(prhs[14]);
    const double normInt[3] = {temp[0], temp[1], temp[2]};
    
    double outMat[3][3] = {};
    for(size_t indNode = 0; indNode < numNodes; ++indNode)
    {
        double tmpMat[3][3] = {};
        const double nodoInt[3] = {intNodes[indNode], intNodes[indNode + numNodes], intNodes[indNode + (2*numNodes)]};
        const double pesoInt = intWeights[indNode];
        const double vettX[3] = {nodeExt[0] - nodoInt[0], nodeExt[1] - nodoInt[1], nodeExt[2] - nodoInt[2]};
        const double lungX = norm2(vettX);
        
        //Aggiunta delle delle 9 componenti del nucleo KR
        nucleoKR(tmpMat, vettX, lungX, currT, cP, cS, lambda, mu, rho, normInt, matCoeff, vettVMS);
        //Aggiunta delle delle 9 componenti dei nuclei KT e KL
        nucleoKLT(tmpMat, vettX, normInt, lungX, currT, nodoInt, matCoeff, vetCoeff, indSM, cP, cS, lambda, mu);
        
        for (size_t r = 0; r < 3; ++r)
            for (size_t c = 0; c < 3; ++c)
                outMat[r][c] += pesoInt * tmpMat[r][c];
    }

    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* outPtr = mxGetPr(plhs[0]);

    for (size_t c = 0; c < 3; ++c)
        for (size_t r = 0; r < 3; ++r)
            outPtr[c*3 + r] = outMat[r][c];
}

// --- FUNCTIONS DEFINITIONS ---
static inline double baseFunctionSM(const double* restrict nodeInt, const double (* restrict matCoeff)[3], const double* restrict vetCoeff, const size_t indSM)
{
    double leftVector[3];
    leftVector[0] = vetCoeff[0] + matCoeff[0][0] * nodeInt[0] + matCoeff[0][1] * nodeInt[1] + matCoeff[0][2] * nodeInt[2];
    leftVector[1] = vetCoeff[1] + matCoeff[1][0] * nodeInt[0] + matCoeff[1][1] * nodeInt[1] + matCoeff[1][2] * nodeInt[2];
    leftVector[2] = vetCoeff[2] + matCoeff[2][0] * nodeInt[0] + matCoeff[2][1] * nodeInt[1] + matCoeff[2][2] * nodeInt[2];
    double rightVector[3] = {0.0, 0.0, 0.0};
    rightVector[indSM - 1] = 1.0;
    return dotProd3D(leftVector, rightVector);
}

static inline void nucleoKLT(double (* restrict nuKLT)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double* restrict nodeInt, const double (* restrict matCoeff)[3], const double* restrict vetCoeff, const size_t indSM, 
               const double cP, const double cS, const double lambda, const double mu)
{
    double nucleoTemp[3][3] = {0};
            
    // Calcolo nucleo KL
    nucleoKL(nucleoTemp, x, n, r, t, cP, cS, lambda, mu);
    // Calcolo nucleo KT
    nucleoKT(nucleoTemp, x, n, r, t, cP, cS, lambda, mu);
    
    // Calcolo funzione di base
    const double baseFunctionValue = baseFunctionSM(nodeInt, matCoeff, vetCoeff, indSM);
    
    // Applicazione componente relativa alla funzione di base
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++) 
            nuKLT[i][j] += baseFunctionValue * nucleoTemp[i][j];
}

static inline void nucleoKL(double (* restrict nuKL)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    const double elst1Factor = lambda/(lambda + mu);
    const double elst2Factor = mu / (lambda + mu);
    const double PWaveFactor = (t - (r/cP)) > 0 ? (t - (r/cP)) : 0.0;
    const double SWaveFactor = (t - (r/cS)) > 0 ? (t - (r/cS)) : 0.0;
    const double wavesFactor = PWaveFactor - SWaveFactor;

    const double invP3Factor = 1 / pow3(r);

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            nuKL[i][j] += (elst1Factor * x[i] * n[j] + elst2Factor * x[j] * n[i]) * invP3Factor * wavesFactor;
}

static inline void nucleoKT(double (* restrict nuKT)[3], const double* restrict x, const double* restrict n, const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    const double elst1Factor = lambda/(lambda + mu);
    const double elst2Factor = mu / (lambda + mu);
    const double PWaveFactor = (t - (r/cP)) > 0 ? (1 / cP) : 0.0;
    const double SWaveFactor = (t - (r/cS)) > 0 ? (1 / cS) : 0.0;
    const double wavesFactor = PWaveFactor - SWaveFactor;

    const double invP2Factor = 1 / pow2(r);

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            nuKT[i][j] += (elst1Factor * x[i] * n[j] + elst2Factor * x[j] * n[i]) * invP2Factor * wavesFactor;
}

static inline void nucleoKR(double (* restrict nuKR)[3], const double* restrict x, const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho,
                                        const double* restrict normInt, const double (* restrict matCoeff)[3], const double* restrict vettVMS)
{
    double nucleoCurr[3][3] = {0};
            
    // Ciclo sull'indice j dei nuclei nuKRj
    for(size_t j = 0; j < 3; j++)
    {   
         // Calcolo nucleo corrente
         nucleoKRj(nucleoCurr, x, r, t, cP, cS, lambda, mu, rho, j);
     
         // Ciclo sulle 9 componenti
         for(size_t i = 0; i < 3; i++)
            for(size_t k = 0; k < 3; k++) 
                nuKR[i][k] -= vettVMS[j] * nucleoCurr[i][k];
    }
}

static inline void nucleoKRj(double (* restrict nuKRj)[3], const double* restrict x, const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const size_t j)
{
    const double elst3Factor = 2 * mu / rho;
    const double elst4Factor = lambda * mu / (rho * (lambda + mu));

    const double PWav1Factor = (t - (r/cP)) > 0 ? (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 : 0.0;
    const double SWav1Factor = (t - (r/cS)) > 0 ? (t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 : 0.0;
    const double wave1Factor = SWav1Factor - PWav1Factor;

    const double PWav2Factor = (t - (r/cP)) > 0 ? (t - (r/cP)) / (cP*cP) : 0.0;
    const double SWav2Factor = (t - (r/cS)) > 0 ? (t - (r/cS)) / (cS*cS) : 0.0;
    const double wave2Factor = SWav2Factor - PWav2Factor;

    const double invP1Factor = 1 / r;
    const double invP3Factor = pow3(invP1Factor);
    const double invP5Factor = pow5(invP1Factor);

    const double e3_i3_w1 = elst3Factor * invP3Factor * wave1Factor;
    const double e3_i3_w2 = elst3Factor * invP3Factor * wave2Factor;
    const double e3_i5_w1 = elst3Factor * 3.0 * invP5Factor * wave1Factor;
    const double e4_i1_w2 = elst4Factor * invP1Factor * wave2Factor;
    
    for (size_t r = 0; r < 3; ++r)
        for (size_t c = 0; c < 3; ++c) 
        {
            double sum_w1 = 0.0;
            double sum_w2 = 0.0;

            for (size_t k = 0; k < 3; ++k) 
            {
                if (k == r) continue;
    
                double x_product = x[k] * x[c];
                sum_w1 += tensorE[j][r][k] * x_product;
                sum_w2 += tensorE[j][r][k] * x_product;
            }
    
            double val = (sum_w1 * e3_i5_w1) + (sum_w2 * e3_i3_w2);

            if(c != r)
                val += (tensorE[j][c][r] * e4_i1_w2) - (tensorE[j][r][c] * e3_i3_w1);

            nuKRj[r][c] = val;
        }
}