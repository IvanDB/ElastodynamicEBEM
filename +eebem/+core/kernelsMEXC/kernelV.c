#include "mex.h"
#include "stdint.h"
#include <math.h>

static inline void nucleo(double (* restrict nu)[3], const double* restrict x, const double r, const double t, const double cP, const double cS);
static inline double norm2(const double v[3]) {return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);}

static inline double pow2(const double x) {return x * x;}
static inline double pow3(const double x) {return x * x * x;}
static inline double pow5(const double x) {const double x2 = pow2(x); return x2 * x2 * x;}

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
    
    double outMat[3][3] = {0.}, tmpMat[3][3] = {0.};

    for(size_t indNode = 0; indNode < numNodes; ++indNode)
    {
        const double nodoInt[3] = {intNodes[indNode], intNodes[indNode + numNodes], intNodes[indNode + (2*numNodes)]};
        const double pesoInt = intWeights[indNode];

        const double vettX[3] = {nodoInt[0] - nodeExt[0], nodoInt[1] - nodeExt[1], nodoInt[2] - nodeExt[2]};
        const double lungX = norm2(vettX);
        
        nucleo(tmpMat, vettX, lungX, currT, cP, cS);
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

static inline void nucleo(double (* restrict nu)[3], const double* restrict x, const double r, const double t, const double cP, const double cS)
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