__device__ void nucleoKL(double nuKL[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
__device__ void nucleoKT(double nuKT[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
__device__ void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, const double cP, const double cS, const double lambda, const double mu);

__device__ void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j);
__device__ void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const double vettRcoeffInt[3]);

__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen);
__device__ double dotProd3D(const double vettA[3], const double vettB[3]);

__global__ void kernelKinternal(double *matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double const4PiDeltaT,
                                  const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                                  const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz, const int numPointSing,
                                  const double *vertsT, const double *areeT, const double *normT, const int offsetZ, const int numBlocks, const double maxLen,
                                  const double *matCoeff, const double *vetCoeff)    
{
    //Evito i blocchi singolari
    if(blockIdx.x == blockIdx.y)
        return;

    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Controllo condizione teorica necessità di calcolo dell'intero blocco
    if (isBlockNull(vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
        return;

    extern __shared__ double matrixSubBlock[][9][9];
    int k, i, j, l, n, m, p;

    const unsigned int sharedBaseInd = threadIdx.x;

    // Inizializzazione shared memory
    for (i = 0; i < 9; i++)
        for (j = 0; j < 9; j++)
            matrixSubBlock[sharedBaseInd][i][j] = 0;

    //Dichiarazione e inizializzazione tensore di permutazione E
    const int tensorE[3][3][3] = {{{0, 0, 0}, {0, 0, 1}, {0, -1, 0}}, {{0, 0, -1}, {0, 0, 0}, {1, 0, 0}}, {{0, 1, 0}, {-1, 0, 0}, {0, 0, 0}}};

    //Dichiarazione variabili
    double vertsTempExt[3][3];
    double vertsTempInt[3][3];
    double nodoTempExt[3];
    double nodoTempInt[3];

    //Lettura vertici triangolo esterno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTempExt[i][j] = vertsT[9*blockIdx.x + 3*j + i];   

    //Estrazione vettore normale corrente
    double normIntCurr[3];
    for (i = 0; i < 3; i++)
        normIntCurr[i] = normT[3*blockIdx.y + i];
    // Lettura vertici triangolo interno
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsTempInt[i][j] = vertsT[9*blockIdx.y + 3*j + i];

    //Ciclo sui nodi della sottoregione corrente
    for (p = 0; p < numPointSing; p++)
    {
        //Calcolo peso nodo GHC corrente
        double pesoInt = stdGHCw[p] * areeT[blockIdx.y];
    
        // Lettura coordinate nodo GHC corrente su triangolo stardard
        nodoTempInt[0] = stdGHCnx[threadIdx.x*numPointSing + p];
        nodoTempInt[1] = stdGHCny[threadIdx.x*numPointSing + p];
        nodoTempInt[2] = stdGHCnz[threadIdx.x*numPointSing + p];
    
        //Mappaggio nodo GHC corrente su triangolo interno corrente
        double nodoGHCcurr[3];
        nodoGHCcurr[0] = nodoTempInt[0] * vertsTempInt[0][0] + nodoTempInt[1] * vertsTempInt[1][0] + nodoTempInt[2] * vertsTempInt[2][0];
        nodoGHCcurr[1] = nodoTempInt[0] * vertsTempInt[0][1] + nodoTempInt[1] * vertsTempInt[1][1] + nodoTempInt[2] * vertsTempInt[2][1];
        nodoGHCcurr[2] = nodoTempInt[0] * vertsTempInt[0][2] + nodoTempInt[1] * vertsTempInt[1][2] + nodoTempInt[2] * vertsTempInt[2][2];
    
    
        //Lettura matrice e vettore dei coefficienti triangolo esterno
        double matCoefCurrExt[3][3];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matCoefCurrExt[i][j] = matCoeff[9*blockIdx.x + 3*j + i];
    
        double vetCoefCurrExt[3];
        for (i = 0; i < 3; i++)
            vetCoefCurrExt[i] = vetCoeff[3*blockIdx.x + i];
    
        //Lettura matrice e vettore dei coefficienti triangolo interno
        double matCoefCurrInt[3][3];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matCoefCurrInt[i][j] = matCoeff[9*blockIdx.y + 3*j + i];
    
        double vetCoefCurrInt[3];
        for (i = 0; i < 3; i++)
            vetCoefCurrInt[i] = vetCoeff[3*blockIdx.y + i];
    
    
        //Ciclo sui numPointExt nodi del triangolo esterno
        for(l = 0; l < numPointExt; l++)
        {
            // Calcolo peso nodo GH triangolo esterno
            double pesoExt = stdGHw[l] * areeT[blockIdx.x];
    
            //Lettura nodo GH corrente su triangolo standard
            nodoTempExt[0] = stdGHnx[l];
            nodoTempExt[1] = stdGHny[l];
            nodoTempExt[2] = stdGHnz[l];
    
            //Mappaggio nodo GH corrente su triangolo esterno
            double nodoGHcurr[3];
            nodoGHcurr[0] = nodoTempExt[0] * vertsTempExt[0][0] + nodoTempExt[1] * vertsTempExt[1][0] + nodoTempExt[2] * vertsTempExt[2][0];
            nodoGHcurr[1] = nodoTempExt[0] * vertsTempExt[0][1] + nodoTempExt[1] * vertsTempExt[1][1] + nodoTempExt[2] * vertsTempExt[2][1];
            nodoGHcurr[2] = nodoTempExt[0] * vertsTempExt[0][2] + nodoTempExt[1] * vertsTempExt[1][2] + nodoTempExt[2] * vertsTempExt[2][2];
    
            //Calcolo coordinate vettore differenza
            double point[3];
            point[0] = nodoGHcurr[0] - nodoGHCcurr[0];
            point[1] = nodoGHcurr[1] - nodoGHCcurr[1];
            point[2] = nodoGHcurr[2] - nodoGHCcurr[2];
    
            //Calcolo norma vettore differenza
            double pointNorm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
    
            //Inizializzazione variabili
            double istTemp = 0;
            double tempValuesLT[3][3];
            double tempValuesR[3][3];
            const int coeffNucleo[4] = {-1, 3, -3, 1};
            const int coeffTemp[4] = {-2, -1, 0, 1};
    
            //Ciclo sui 4 istanti temporali
            for (k = 0; k < 4; k++)
            {
                // Calcolo istante temporale corrente
                istTemp = deltaT * (offsetZ + double(blockIdx.z) + coeffTemp[k]);
    
                //Check necessità di calcolo
                if(istTemp <= 0)
                   continue;
    
                for (n = 0; n < 3; n++)
                    for (m = 0; m < 3; m++)
                    {
                        double baseFuncValueExt = matCoefCurrExt[n][0] * nodoGHcurr[0] +  matCoefCurrExt[n][1] * nodoGHcurr[1] + matCoefCurrExt[n][2] * nodoGHcurr[2] + vetCoefCurrExt[n];
                        double baseFuncValueInt = matCoefCurrInt[m][0] * nodoGHCcurr[0] +  matCoefCurrInt[m][1] * nodoGHCcurr[1] + matCoefCurrInt[m][2] * nodoGHCcurr[2] + vetCoefCurrInt[m];
                        
                        double vettRcoeffInt[3] = {+ normIntCurr[1] * matCoefCurrInt[m][2] - normIntCurr[2] * matCoefCurrInt[m][1],
                                                   - normIntCurr[0] * matCoefCurrInt[m][2] + normIntCurr[2] * matCoefCurrInt[m][0],
                                                   + normIntCurr[0] * matCoefCurrInt[m][1] - normIntCurr[1] * matCoefCurrInt[m][0]};
    
                        
    
                        //Inizializzazione componenti nucleo totale
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 3; j++)
                            {
                                tempValuesLT[i][j] = 0;
                                tempValuesR[i][j] = 0;
                            }
            
                        //Aggiunta delle delle 9 componenti del nucleo KR
                        nucleoKR(tempValuesR, point, pointNorm, istTemp, velP, velS, lambda, mu, rho, tensorE, vettRcoeffInt);
                        
                        //Aggiunta delle delle 9 componenti dei nuclei KT e KL
                        nucleoKT(tempValuesLT, point, normIntCurr, pointNorm, istTemp, velP, velS, lambda, mu);
                        nucleoKL(tempValuesLT, point, normIntCurr, pointNorm, istTemp, velP, velS, lambda, mu);
                        
                        //Somma pesata dei valori del nucleo alla shared memory
                        for (i = 0; i < 3; i++)
                            for (j = 0; j < 3; j++)
                                matrixSubBlock[sharedBaseInd][3*n + i][3*m + j] += pesoExt * pesoInt * coeffNucleo[k] * (baseFuncValueExt * (tempValuesR[i][j] + (baseFuncValueInt * tempValuesLT[i][j]))); 
                    }
            }
        }
    }

    //Sync prima di inziare la riduzione
    __syncthreads();

    unsigned int xDim, sharedOffInd;

    //Iterazione per scendere ad 1 in x (x è sempre potenza di 2)
    xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = threadIdx.x + xDim;
        if(threadIdx.x < xDim)
            for (i = 0; i < 9; i++)
                for (j = 0; j < 9; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    // Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0)
        for (i = 0; i < 9; i++)
            for (j = 0; j < 9; j++)
            {
                ind = 81*gridDim.x*gridDim.y*blockIdx.z + 9*gridDim.x*(9*blockIdx.y + j) + 9*blockIdx.x + i;
                //matrix[9*blockIdx.x + i][9*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiDeltaT;
            }

    
}

__device__ bool isBlockNull(const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
{
    int i, j;

    double vertsS[3][3];
    double vertsF[3][3];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];
            vertsF[i][j] = vertsT[9*blockIdx.y + 3*j + i];
        }
        

    double distMax = 0;
    double distMin = 100; //Mettere valore vicino a +Inf;
    double distCurr;

    double vettDist[3];

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            vettDist[0] = vertsS[i][0] - vertsF[j][0];
            vettDist[1] = vertsS[i][1] - vertsF[j][1];
            vettDist[2] = vertsS[i][2] - vertsF[j][2];
            distCurr = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]);
            if (distCurr > distMax)
                distMax = distCurr;
            if (distCurr < distMin)
                distMin = distCurr;
        }
    
    return ((indTemp - 2) * cS * deltaT > distMax) || ((indTemp + 1) * cP * deltaT < distMin - maxLen);
}


__device__ double dotProd3D(const double vettA[3], const double vettB[3])
{
    return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];
}


__device__ void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
                                    const double cP, const double cS, const double lambda, const double mu)
{
    // Calcolo nucleo KL
    nucleoKL(nuKLT, x, n, r, t, cP, cS, lambda, mu);

    // Calcolo nuckeo KT
    nucleoKT(nuKLT, x, n, r, t, cP, cS, lambda, mu);
}


__device__ void nucleoKL(double nuKL[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    nuKL[0][0] += (lambda/(lambda + mu)) * x[0] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                                + dotProd3D(x, n) / pow(r, 3) * (t - (r/cS)) * ((t - (r/cS)) > 0);

    nuKL[0][1] += (lambda/(lambda + mu)) * x[0] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[0][2] += (lambda/(lambda + mu)) * x[0] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[1][0] += (lambda/(lambda + mu)) * x[1] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[1][1] += (lambda/(lambda + mu)) * x[1] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                                + dotProd3D(x, n) / pow(r, 3) * (t - (r/cS)) * ((t - (r/cS)) > 0);

    nuKL[1][2] += (lambda/(lambda + mu)) * x[1] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][0] += (lambda/(lambda + mu)) * x[2] * n[0] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[0] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][1] += (lambda/(lambda + mu)) * x[2] * n[1] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[1] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0));

    nuKL[2][2] += (lambda/(lambda + mu)) * x[2] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                        +  (mu / (lambda + mu)) * x[2] * n[2] / pow(r, 3) * ((t - (r/cP)) * ((t - (r/cP)) > 0) - (t - (r/cS)) * ((t - (r/cS)) > 0))
                                + dotProd3D(x, n) / pow(r, 3) * (t - (r/cS)) * ((t - (r/cS)) > 0);
}

__device__ void nucleoKT(double nuKT[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu)
{
    nuKT[0][0] += (lambda/(lambda + mu)) * x[0] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                                + dotProd3D(x, n) / pow(r, 2) * (((t - (r/cS)) > 0) / cS);

    nuKT[0][1] += (lambda/(lambda + mu)) * x[0] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[0][2] += (lambda/(lambda + mu)) * x[0] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[1][0] += (lambda/(lambda + mu)) * x[1] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[1][1] += (lambda/(lambda + mu)) * x[1] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                                + dotProd3D(x, n) / pow(r, 2) * (((t - (r/cS)) > 0) / cS);

    nuKT[1][2] += (lambda/(lambda + mu)) * x[1] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][0] += (lambda/(lambda + mu)) * x[2] * n[0] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[0] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][1] += (lambda/(lambda + mu)) * x[2] * n[1] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[1] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS));

    nuKT[2][2] += (lambda/(lambda + mu)) * x[2] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                        +  (mu / (lambda + mu)) * x[2] * n[2] / pow(r, 2) * ((((t - (r/cP)) > 0) / cP) - (((t - (r/cS)) > 0) / cS))
                                + dotProd3D(x, n) / pow(r, 2) * (((t - (r/cS)) > 0) / cS);
}

__device__ void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const double vettRcoeffInt[3])
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
                nuKR[i][k] -= nucleoCurr[i][k] * vettRcoeffInt[j];
    }
}

__device__ void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j)
{
    nuKRj[0][0] = (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[0] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[0] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP)))
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][0][2] * (3 * x[2] * x[1] / pow(r, 5)))
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[1] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[0][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][0] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (3 * x[1] * x[2] / pow(r, 5)) + tensorE[j][0][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][0][1] * (x[1] * x[2] / pow(r, 3)) + tensorE[j][0][2] * (x[2] * x[2] / pow(r, 3)))
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][1][2] * (3 * x[2] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][1] = (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[1] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[1][2] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][2][1] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][1][2] * (3 * x[2] * x[2] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][1][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][1][2] * (x[2] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][0] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][0][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[0] / pow(r, 5) - 1 / pow(r, 3)) + tensorE[j][2][1] * (3 * x[1] * x[0] / pow(r, 5)))
                                                            * ((t - (r/cS)) * (t - (r/cS)) * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0))
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[0] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[0] / pow(r, 3)))
                                                            * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][1] = (lambda * mu / (rho * (lambda+mu))) * tensorE[j][1][2] / r * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((r - (r/cP)) > 0) / (cP*cP))) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[1] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[1] / pow(r, 5) - 1 / pow(r, 3))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[1] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[1] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));

    nuKRj[2][2] = (2 * mu / rho) * (tensorE[j][2][0] * (3 * x[0] * x[2] / pow(r, 5)) + tensorE[j][2][1] * (3 * x[1] * x[2] / pow(r, 5))) 
                                                           * ((t - (r/cS)) * (t - (r/cS))  * (t + (2*r/cS)) / 6 * ((t - (r/cS)) > 0) - (t - (r/cP)) * (t - (r/cP)) * (t + (2*r/cP)) / 6 * ((t - (r/cP)) > 0)) 
                            + (2 * mu / rho) * (tensorE[j][2][0] * (x[0] * x[2] / pow(r, 3)) + tensorE[j][2][1] * (x[1] * x[2] / pow(r, 3))) 
                                                           * (((t - (r/cS)) * ((t - (r/cS)) > 0) / (cS*cS)) - ((t - (r/cP)) * ((t - (r/cP)) > 0) / (cP*cP)));
}