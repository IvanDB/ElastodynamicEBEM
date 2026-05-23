#include "mex.h"
#include "stdint.h"

void nucleoKL(double nuKL[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
void nucleoKT(double nuKT[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, const double cP, const double cS, const double lambda, const double mu);

void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j);
void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const double normInt[3], const double matCoeff[3][3], const double vettVMS[3]);

double dotProd3D(const double vettA[3], const double vettB[3]) {return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];}
double norm2(const double v[3]) {return sqrt(dotProd3D(v, v));}

double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs == 0)
        return;

    double* temp = 0x0;
    
    temp = mxGetPr(prhs[0]);
    int numNodes = (int) temp[0];

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

    temp = mxGetPr(prhs[7]);
    double lambda = temp[0];
    
    temp = mxGetPr(prhs[8]);
    double mu = temp[0];
    
    temp = mxGetPr(prhs[9]);
    double rho = temp[0];

    temp = mxGetPr(prhs[10]);
    double vettVMS[3] = {temp[0], temp[1], temp[2]};

    temp = mxGetPr(prhs[11]);
    double matCoeff[3][3]  = {{temp[0], temp[3], temp[6]}, 
                              {temp[1], temp[4], temp[7]}, 
                              {temp[2], temp[5], temp[8]}};

    temp = mxGetPr(prhs[12]);
    double vetCoeff[3] = {temp[0], temp[1], temp[2]};

    temp = mxGetPr(prhs[13]);
    int indSM = (int) temp[0];

    temp = mxGetPr(prhs[14]);
    double normInt[3] = {temp[0], temp[1], temp[2]};

    //Dichiarazione e inizializzazione tensore di permutazione E
    const int tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, 
                                  {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, 
                                  {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

    
    double outMat[3][3] = {};

    for(int indNode = 0; indNode < numNodes; ++indNode)
    {
        double tmpMat[3][3] = {};
        double nodoInt[3] = {intNodes[indNode], intNodes[indNode + numNodes], intNodes[indNode + (2*numNodes)]};
        double pesoInt = intWeights[indNode];

        double vettX[3] = {nodeExt[0] - nodoInt[0], nodeExt[1] - nodoInt[1], nodeExt[2] - nodoInt[2]};
        double lungX = norm2(vettX);
        
        //Aggiunta delle delle 9 componenti del nucleo KR
        nucleoKR(tmpMat, vettX, lungX, currT, cP, cS, lambda, mu, rho, tensorE, normInt, matCoeff, vettVMS);

        //Aggiunta delle delle 9 componenti dei nuclei KT e KL
        nucleoKLT(tmpMat, vettX, normInt, lungX, currT, nodoInt, matCoeff, vetCoeff, indSM, cP, cS, lambda, mu);

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


double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM)
{
    double leftVector[3];
    leftVector[0] = vetCoeff[0] + matCoeff[0][0] * nodeInt[0] + matCoeff[0][1] * nodeInt[1] + matCoeff[0][2] * nodeInt[2];
    leftVector[1] = vetCoeff[1] + matCoeff[1][0] * nodeInt[0] + matCoeff[1][1] * nodeInt[1] + matCoeff[1][2] * nodeInt[2];
    leftVector[2] = vetCoeff[2] + matCoeff[2][0] * nodeInt[0] + matCoeff[2][1] * nodeInt[1] + matCoeff[2][2] * nodeInt[2];

    double rightVector[3] = {0.0, 0.0, 0.0};
    rightVector[indSM - 1] = 1.0;
    return dotProd3D(leftVector, rightVector);
}

void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
               const double cP, const double cS, const double lambda, const double mu)
{
    int i, j;
    double nucleoTemp[3][3];
    //Inizializzazione nucleo temporaneo
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) 
            nucleoTemp[i][j] = 0;

    // Calcolo nucleo KL
    nucleoKL(nucleoTemp, x, n, r, t, cP, cS, lambda, mu);

    // Calcolo nuckeo KT
    nucleoKT(nucleoTemp, x, n, r, t, cP, cS, lambda, mu);
    
    // Calcolo funzione di base
    const double baseFunctionValue = baseFunctionSM(nodeInt, matCoeff, vetCoeff, indSM);

    //Applicazione componente relativa alla funzione di base
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) 
            nuKLT[i][j] += baseFunctionValue * nucleoTemp[i][j];
}


void nucleoKL(double nuKL[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    nuKL[0][0] += (lambda/(lambda + mu)) * x[0] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[0][1] += (lambda/(lambda + mu)) * x[0] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[0][2] += (lambda/(lambda + mu)) * x[0] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[1][0] += (lambda/(lambda + mu)) * x[1] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[1][1] += (lambda/(lambda + mu)) * x[1] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[1][2] += (lambda/(lambda + mu)) * x[1] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][0] += (lambda/(lambda + mu)) * x[2] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][1] += (lambda/(lambda + mu)) * x[2] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][2] += (lambda/(lambda + mu)) * x[2] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));
}

void nucleoKT(double nuKT[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    nuKT[0][0] += (lambda/(lambda + mu)) * x[0] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[0][1] += (lambda/(lambda + mu)) * x[0] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[0][2] += (lambda/(lambda + mu)) * x[0] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[1][0] += (lambda/(lambda + mu)) * x[1] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[1][1] += (lambda/(lambda + mu)) * x[1] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[1][2] += (lambda/(lambda + mu)) * x[1] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][0] += (lambda/(lambda + mu)) * x[2] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][1] += (lambda/(lambda + mu)) * x[2] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][2] += (lambda/(lambda + mu)) * x[2] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));
}


void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3],
                                        const double normInt[3], const double matCoeff[3][3], const double vettVMS[3])
{
    int j, i, k;
    double nucleoCurr[3][3];
    // Inizializzazione nucleo temporaneo
    for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++) 
            nucleoCurr[i][k] = 0;

    // Ciclo sull'indice j dei nuclei nuKRj
    for(j = 0; j < 3; j++)
    {   
         // Calcolo nucleo corrente
         nucleoKRj(nucleoCurr, x, r, t, cP, cS, lambda, mu, rho, tensorE, j);
     
         // Ciclo sulle 9 componenti
         for (i = 0; i < 3; i++)
             for (k = 0; k < 3; k++) 
                 nuKR[i][k] -= vettVMS[j] * nucleoCurr[i][k];
    }
}

void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j)
{
    nuKRj[0][0] = (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[0] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[0] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][0][2] * (3 * x[2] * x[1] / pow(r, 5)))
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[1] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[2] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[2] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[2] / pow(r, 3)))
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][1][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][1] = (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[1] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][2][1] * (3 * x[1] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][2] = (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[2] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));
}

