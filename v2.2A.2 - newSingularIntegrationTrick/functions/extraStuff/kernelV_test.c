#include "mex.h"
#include "stdint.h"

void nucleo(double out[3][3], const double x[3], const double r, const double t, const double cP, const double cS);
double norm2(const double v[3]) {return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs == 0)
        return;

    double* temp = 0x0;
    
    temp = mxGetPr(prhs[0]);
    uint16_t numNodes = (uint16_t) temp[0];

    temp = mxGetPr(prhs[1]);
    double* intNodes = temp;

    temp = mxGetPr(prhs[2]);
    double* intWeights = temp;
    
    temp = mxGetPr(prhs[3]);
    double nodeExt[3] = {temp[0], temp[1], temp[2]};

    temp = mxGetPr(prhs[4]);
    double currT = temp[0];

    temp = mxGetPr(prhs[5]);
    double cP = temp[0];

    temp = mxGetPr(prhs[6]);
    double cS = temp[0];
    
    double outMat[3][3] = {}, tmpMat[3][3] = {};

    for(uint16_t indNode = 0; indNode < numNodes; ++indNode)
    {
        double nodoInt[3] = {intNodes[indNode], intNodes[indNode + numNodes], intNodes[indNode + (2*numNodes)]};
        double pesoInt = intWeights[indNode];

        double vettX[3] = {nodoInt[0] - nodeExt[0], nodoInt[1] - nodeExt[1], nodoInt[2] - nodeExt[2]};
        double lungX = norm2(vettX);
        
        nucleo(tmpMat, vettX, lungX, currT, cP, cS);
        for (uint8_t c = 0; c < 3; ++c)
            for (uint8_t r = 0; r < 3; ++r)
                outMat[r][c] += pesoInt * tmpMat[r][c];
    }

    plhs[0] = mxCreateDoubleMatrix(3, 3, mxREAL);
    double* outPtr = mxGetPr(plhs[0]);

    for (uint8_t c = 0; c < 3; ++c)
        for (uint8_t r = 0; r < 3; ++r)
            outPtr[c*3 + r] = outMat[r][c];
    
}

void nucleo(double nu[3][3], const double x[3], const double r, const double t, const double cP, const double cS)
{
    nu[0][0] = ((x[0] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[0] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[0] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[0][1] = ((x[0] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[0] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[0] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[0][2] = ((x[0] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[0] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[0] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[1][0] = ((x[1] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[1] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[1] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[1][1] = ((x[1] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[1] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[1] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[1][2] = ((x[1] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[1] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[1] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][0] = ((x[2] * x[0]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[2] * x[0]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[2] * x[0]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][1] = ((x[2] * x[1]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ( - ((x[2] * x[1]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * (- 3 * ((x[2] * x[1]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));

    nu[2][2] = ((x[2] * x[2]) / (pow(r, 3))) * (1 / (cP*cP)) * (t - (r/cP) > 0)
                    +  ((1/r) - ((x[2] * x[2]) / (pow(r, 3)))) * (1 / (cS*cS)) * (t - (r/cS) > 0)
                    - 0.5 * ((1/(pow(r, 3))) - 3 * ((x[2] * x[2]) / (pow(r, 5)))) *
                                        ((cP*cP * t*t - r*r) / (cP*cP) * (t - (r/cP) > 0)
                                          - (cS*cS * t*t - r*r) / (cS*cS) * (t - (r/cS) > 0));
}