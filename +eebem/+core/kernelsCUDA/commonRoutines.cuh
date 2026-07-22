#include <float.h>

//Macro
#ifdef __CUDACC__
    #define __eebemDevice __device__
#else
    #define __eebemDevice
#endif

//Basic math functions
static __eebemDevice inline double pow2(const double x) 
{
    return x * x;
}

static __eebemDevice inline double pow3(const double x) 
{
    return x * x * x;
}

static __eebemDevice inline double pow5(const double x) 
{
    const double x2 = pow2(x);
    return x2 * x2 * x;
}

static __eebemDevice inline double norm2(const double v[3]) 
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static __eebemDevice inline double dotProd3D(const double v1[3], const double v2[3])
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

//Permutation tensor
static __eebemDevice const double tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, 
                                                      {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}},
                                                      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

//Auxiliary functions
static __eebemDevice inline double baseFunctionSM(const double* __restrict__ nodeInt, const double (*__restrict__ matCoeff)[3], const double* __restrict__ vetCoeff, const int indSM)
{
    const double leftVector[3] = {vetCoeff[0] + matCoeff[0][0] * nodeInt[0] + matCoeff[0][1] * nodeInt[1] + matCoeff[0][2] * nodeInt[2],
                                  vetCoeff[1] + matCoeff[1][0] * nodeInt[0] + matCoeff[1][1] * nodeInt[1] + matCoeff[1][2] * nodeInt[2],
                                  vetCoeff[2] + matCoeff[2][0] * nodeInt[0] + matCoeff[2][1] * nodeInt[1] + matCoeff[2][2] * nodeInt[2]};

    double rightVector[3] = {0.};
    rightVector[indSM - 1] = 1.0;
    return dotProd3D(leftVector, rightVector);
}


//Sub-kernels
static __eebemDevice inline void nuV(double (* __restrict__ nu)[3], const double* __restrict__ x, const double r, const double t, const double cP, const double cS)
{
    const double PWav1Factor = (t - (r/cP)) > 0 ? 1 / pow2(cP) : 0.0;
    const double SWav1Factor = (t - (r/cS)) > 0 ? 1 / pow2(cS) : 0.0;
    const double wave1Factor = PWav1Factor - SWav1Factor;

    const double PWav2Factor = (t - (r/cP)) > 0 ? pow2(t) - pow2(r/cP) : 0.0;
    const double SWav2Factor = (t - (r/cS)) > 0 ? pow2(t) - pow2(r/cS) : 0.0;
    const double wave2Factor = 0.5 * (PWav2Factor - SWav2Factor);

    const double invP1Factor = 1 / r;
    const double invP3Factor = pow3(invP1Factor);
    const double invP5Factor = pow5(invP1Factor);
    
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
        {
            const double xi_xj_P3 = x[i] * x[j] * invP3Factor;
            const double xi_xj_P5 = x[i] * x[j] * invP5Factor;

            nu[i][j] = xi_xj_P3 * wave1Factor 
                + 3. * xi_xj_P5 * wave2Factor;

            if(i == j)
                nu[i][j] += invP1Factor * SWav1Factor - invP3Factor * wave2Factor;
        }
}


static __eebemDevice inline void nuKL(double (* __restrict__ nuKL)[3], const double* __restrict__ x, const double* __restrict__ n, const double r, const double t, 
                                      const double cP, const double cS, const double lambda, const double mu)
{
    const double elst1Factor = lambda/(lambda + mu);
    const double elst2Factor = mu / (lambda + mu);
    const double PWaveFactor = (t - (r/cP)) > 0 ? (t - (r/cP)) : 0.0;
    const double SWaveFactor = (t - (r/cS)) > 0 ? (t - (r/cS)) : 0.0;
    const double wavesFactor = PWaveFactor - SWaveFactor;

    const double invP3Factor = 1 / pow3(r);
    const double diagAddTerm = dotProd3D(x, n) * invP3Factor * SWaveFactor;

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            nuKL[i][j] += (elst1Factor * x[i] * n[j] + elst2Factor * x[j] * n[i]) * invP3Factor * wavesFactor;
            if(i == j)
                nuKL[i][j] += diagAddTerm;
        }
}

static __eebemDevice inline void nuKT(double (* __restrict__ nuKT)[3], const double* __restrict__ x, const double* __restrict__ n, const double r, const double t, 
                                      const double cP, const double cS, const double lambda, const double mu)
{
    const double elst1Factor = lambda/(lambda + mu);
    const double elst2Factor = mu / (lambda + mu);
    const double PWaveFactor = (t - (r/cP)) > 0 ? (1 / cP) : 0.0;
    const double SWaveFactor = (t - (r/cS)) > 0 ? (1 / cS) : 0.0;
    const double wavesFactor = PWaveFactor - SWaveFactor;

    const double invP2Factor = 1 / pow2(r);
    const double diagAddTerm = dotProd3D(x, n) * invP2Factor * SWaveFactor;

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            nuKT[i][j] += (elst1Factor * x[i] * n[j] + elst2Factor * x[j] * n[i]) * invP2Factor * wavesFactor;
            if(i == j)
                nuKT[i][j] += diagAddTerm;
        }
}

static __eebemDevice inline void nuKLT(double (* __restrict__ nuKLT)[3], const double* __restrict__ x, const double* __restrict__ n, const double r, const double t, 
                                       const double* __restrict__ nodeInt, const double (*__restrict__ matCoeff)[3], const double* __restrict__ vetCoeff, const int indSM, 
                                       const double cP, const double cS, const double lambda, const double mu)
{
    double nuTemp[3][3] = {0.};

    //Add KL kernel
    nuKL(nuTemp, x, n, r, t, cP, cS, lambda, mu);

    //Add KT kernel
    nuKT(nuTemp, x, n, r, t, cP, cS, lambda, mu);
    
    //Base function
    const double baseFunctionValue = baseFunctionSM(nodeInt, matCoeff, vetCoeff, indSM);
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j) 
            nuKLT[i][j] += nuTemp[i][j] * baseFunctionValue;
}

static __eebemDevice inline void nuKRj(double (* __restrict__ nuKRj)[3], const double x[3], const double r, const double t, const double cP, const double cS, 
                                       const double lambda, const double mu, const double rho, const size_t j)
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

    for(size_t r = 0; r < 3; ++r)
        for(size_t c = 0; c < 3; ++c) 
        {
            double sum_w1 = 0.0;
            double sum_w2 = 0.0;

            for(size_t k = 0; k < 3; ++k) 
            {
                if (k == r) continue;

                const double x_product = x[k] * x[c];
                sum_w1 += tensorE[j][r][k] * x_product;
                sum_w2 += tensorE[j][r][k] * x_product;
            }

            double val = (sum_w1 * e3_i5_w1) + (sum_w2 * e3_i3_w2);

            if(c != r)
                val += (tensorE[j][c][r] * e4_i1_w2) - (tensorE[j][r][c] * e3_i3_w1);

            nuKRj[r][c] = val;
        }
}

static __eebemDevice inline void nuKR_linear(double (* __restrict__ nuKR)[3], const double* __restrict__ x, const double r, const double t, 
                                      const double cP, const double cS, const double lambda, const double mu, const double rho,
                                      const double* __restrict__ normInt, const double (* __restrict__ matCoeff)[3], const int indSM)
{
    double nuTemp[3][3] = {0.};

    //Vector Vsm 
    const int zeroBasedIndSM = indSM - 1;
    const double vettRSM[3] = {+ normInt[1] * matCoeff[zeroBasedIndSM][2] - normInt[2] * matCoeff[zeroBasedIndSM][1],
                               - normInt[0] * matCoeff[zeroBasedIndSM][2] + normInt[2] * matCoeff[zeroBasedIndSM][0],
                               + normInt[0] * matCoeff[zeroBasedIndSM][1] - normInt[1] * matCoeff[zeroBasedIndSM][0]};

    //Loop over j index of nuKRj
    for(size_t j = 0; j < 3; ++j)
    {   
        nuKRj(nuTemp, x, r, t, cP, cS, lambda, mu, rho, j);
        
        for(size_t i = 0; i < 3; ++i)
            for(size_t k = 0; k < 3; ++k) 
                nuKR[i][k] -= vettRSM[j] * nuTemp[i][k];
    }
}

static __device__ inline void nuKR_costant(double (* __restrict__ nuKR)[3], const double* __restrict__ tau, const double* __restrict__ x, const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho)
{
    double nuTemp[3][3] = {0.};

    //Loop over j index of nuKRj
    for(size_t j = 0; j < 3; ++j)
    {
        nuKRj(nuTemp, x, r, t, cP, cS, lambda, mu, rho, j);

        for(size_t i = 0; i < 3; ++i)
            for(size_t k = 0; k < 3; ++k)
                nuKR[i][k] += nuTemp[i][k] * tau[j];
    }
}

//Null block check functions
static __eebemDevice inline bool isVBlockNull(const double* __restrict__ vertsT, const double deltaT, const double cP, const double cS, const int indTemp)
{
    double vertsS[3][3];
    double vertsF[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];
            vertsF[i][j] = vertsT[9*blockIdx.y + 3*j + i];
        }

    double distMax = 0;

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            const double vettDist[] = {vertsS[i][0] - vertsF[j][0],
                                       vertsS[i][1] - vertsF[j][1],
                                       vertsS[i][2] - vertsF[j][2]};

            const double distCurr = norm2(vettDist);
            if(distCurr > distMax)
                distMax = distCurr;
        }
    
    return ((indTemp - 1) * cS * deltaT > distMax);
}

static __eebemDevice inline bool isKBlockNull_costant(const double* __restrict__ vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
{
    double vertsS[3][3];
    double vertsF[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];
            vertsF[i][j] = vertsT[9*blockIdx.y + 3*j + i];
        }
        
    double distMax = 0;
    double distMin = DBL_MAX;

    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
        {
            const double vettDist[] = {vertsS[i][0] - vertsF[j][0],
                                       vertsS[i][1] - vertsF[j][1],
                                       vertsS[i][2] - vertsF[j][2]};

            const double distCurr = norm2(vettDist);
            if(distCurr > distMax)
                distMax = distCurr;
            if(distCurr < distMin)
                distMin = distCurr;
        }
    
    return ((indTemp - 2) * cS * deltaT > distMax) || ((indTemp + 1) * cP * deltaT < distMin - maxLen);
}

static __eebemDevice inline bool isKBlockNull_linear(const double* __restrict__ nodesMesh, const double* __restrict__ vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
{
    double vertsS[3][3];
    for(size_t i = 0; i < 3; ++i)
        for(size_t j = 0; j < 3; ++j)
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    double pointF[3]; 
    for(size_t i = 0; i < 3; ++i)
        pointF[i] = nodesMesh[3*blockIdx.y + i];
        
    double distMax = 0;
    double distMin = DBL_MAX;

    for(size_t i = 0; i < 3; ++i)
    {
        const double vettDist[] = {vertsS[i][0] - pointF[0],
                                   vertsS[i][1] - pointF[1],
                                   vertsS[i][2] - pointF[2]};

        const double distCurr = norm2(vettDist);
        if(distCurr > distMax)
            distMax = distCurr;
        if(distCurr < distMin)
            distMin = distCurr;
    }

    return ((indTemp - 2) * cS * deltaT > distMax + maxLen) || ((indTemp + 1) * cP * deltaT < distMin - 2*maxLen);
}