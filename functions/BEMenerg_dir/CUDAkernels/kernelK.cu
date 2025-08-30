__device__ void nucleoKL(double nuKL[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
__device__ void nucleoKT(double nuKT[3][3], const double x[3], const double n[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu);
__device__ void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, const double cP, const double cS, const double lambda, const double mu);

__device__ void nucleoKRj(double nuKRj[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const int j);
__device__ void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3], const double normInt[3], const double matCoeff[3][3], const int indSM);

__device__ bool isBlockNull(const double *nodesMesh, const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen);
__device__ double dotProd3D(const double vettA[3], const double vettB[3]);
__device__ double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM);
__device__ void vectorRSM(double vectorRSMValue[3], const double normInt[3], const double matCoeff[3][3], const int indSM);


__global__ void kernelK(double *matrix, const double deltaT, const double velP, const double velS, const double lambda, const double mu, const double rho, const double const4PiDeltaT,
                          const double *stdGHw, const double *stdGHnx, const double *stdGHny, const double *stdGHnz, const int numPointExt,
                          const double *stdGHCw, const double *stdGHCnx, const double *stdGHCny, const double *stdGHCnz,
                          const double *vertsT, const double *areeT, const double *normT, const int *indSMmatrix, const double *matCoeff, const double *vetCoeff, 
                          const int offsetZ, const int numBlocks, const double *nodesMesh, const double maxLen)  
{
    // //Evito i blocchi diagonali
    // if(blockIdx.x == blockIdx.y)
    //     return;

    //Controllo di non aver sforato l'indice temporale
    if(offsetZ + blockIdx.z >= numBlocks)
        return;

    //Controllo condizione teorica necessità di calcolo dell'intero blocco
    if (isBlockNull(nodesMesh, vertsT, deltaT, velP, velS, offsetZ + blockIdx.z, maxLen))
        return;

    extern __shared__ double matrixSubBlock[][3][3];
    int k, i, j, l, m;

    const unsigned int sharedBaseInd = threadIdx.x * blockDim.y + threadIdx.y;

    // Inizializzazione shared memory
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
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

    //Inizializzazione indice matrice delle flag
    const int baseIndFlag = blockIdx.y * gridDim.x;

    //Ciclo sui triangoli di campo
    for(m = 0; m < gridDim.x; m++)
    { 
        //Salto i casi singolari, si fanno su CPU
        if(m == blockIdx.x)
            continue;

        //Estrazione indice vertice-triangolo corrente
        int indSMcurr = indSMmatrix[baseIndFlag + m];

        //Check indice associato al blocco
        if(indSMcurr == 0)
            continue;

        //Estrazione vettore normale corrente
        double normIntCurr[3];
        for (i = 0; i < 3; i++)
            normIntCurr[i] = normT[3*m + i];

        //Estrazione matrice e vettore dei coefficienti correnti
        double matCoeffCurr[3][3];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
               matCoeffCurr[i][j] = matCoeff[9*m + 3*j + i];

        double vetCoeffCurr[3];
        for (i = 0; i < 3; i++)
            vetCoeffCurr[i] = vetCoeff[3*m + i];

        //Calcolo peso nodo GHC corrente
        double pesoInt = stdGHCw[threadIdx.y] * areeT[m];

        // Lettura vertici triangolo interno
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                vertsTempInt[i][j] = vertsT[9*m + 3*j + i];

        // Lettura coordinate nodo GHC corrente su triangolo stardard
        nodoTempInt[0] = stdGHCnx[threadIdx.x*blockDim.y + threadIdx.y];
        nodoTempInt[1] = stdGHCny[threadIdx.x*blockDim.y + threadIdx.y];
        nodoTempInt[2] = stdGHCnz[threadIdx.x*blockDim.y + threadIdx.y];

        //Mappaggio nodo GHC corrente su triangolo interno corrente
        double nodoGHCcurr[3];
        nodoGHCcurr[0] = nodoTempInt[0] * vertsTempInt[0][0] + nodoTempInt[1] * vertsTempInt[1][0] + nodoTempInt[2] * vertsTempInt[2][0];
        nodoGHCcurr[1] = nodoTempInt[0] * vertsTempInt[0][1] + nodoTempInt[1] * vertsTempInt[1][1] + nodoTempInt[2] * vertsTempInt[2][1];
        nodoGHCcurr[2] = nodoTempInt[0] * vertsTempInt[0][2] + nodoTempInt[1] * vertsTempInt[1][2] + nodoTempInt[2] * vertsTempInt[2][2];

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
            double tempValues[3][3];
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

                //Inizializzazione componenti nucleo totale
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        tempValues[i][j] = 0;

                //Aggiunta delle delle 9 componenti del nucleo KR
                nucleoKR(tempValues, point, pointNorm, istTemp, velP, velS, lambda, mu, rho, tensorE, normIntCurr, matCoeffCurr, indSMcurr);

                //Aggiunta delle delle 9 componenti dei nuclei KT e KL
                nucleoKLT(tempValues, point, normIntCurr, pointNorm, istTemp, nodoGHCcurr, matCoeffCurr, vetCoeffCurr, indSMcurr, velP, velS, lambda, mu);
                
                //Somma pesata dei valori del nucleo alla shared memory
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        matrixSubBlock[sharedBaseInd][i][j] += pesoExt * pesoInt * coeffNucleo[k] * tempValues[i][j]; //pesoExt * pesoInt * coeffNucleo[k] * tempValues[i][j] 
            }
        }
    }

    //Sync prima di inziare la riduzione
    __syncthreads();

    unsigned int xDim, yDim, sharedOffInd;

    //Primo step per scendere a potenza di 2 in y
    yDim = pow(2.0, (int) floor(log2((float) blockDim.y)));
    sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;

    if(threadIdx.y + yDim < blockDim.y)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

    __syncthreads();

    //Iterazione per scendere ad 1 in y
    yDim /= 2;
    while(yDim > 0)
    {
        sharedOffInd = threadIdx.x * blockDim.y + threadIdx.y + yDim;
        if(threadIdx.y < yDim)
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        yDim /= 2;
    }

    //Iterazione per scendere ad 1 in x (x è sempre potenza di 2)
    xDim = blockDim.x/2;
    while(xDim > 0)
    {
        sharedOffInd = (threadIdx.x + xDim) * blockDim.y + threadIdx.y;
        if(threadIdx.x < xDim)
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    matrixSubBlock[sharedBaseInd][i][j] += matrixSubBlock[sharedOffInd][i][j];

        __syncthreads();
        xDim /= 2;
    }

    // Salvataggio dati in memoria globale
    unsigned long ind;
    if(threadIdx.x == 0 && threadIdx.y == 0)
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
            {
                ind = 9*gridDim.x*gridDim.y*blockIdx.z + 3*gridDim.x*(3*blockIdx.y + j) + 3*blockIdx.x + i;
                //matrix[3*blockIdx.x + i][3*blockIdx.y + j][blockIdx.z]
                if(abs(matrixSubBlock[0][i][j]) > pow(10.0, -14))
                    matrix[ind] += matrixSubBlock[0][i][j] / const4PiDeltaT;
            }

    
}

__device__ bool isBlockNull(const double *nodesMesh, const double *vertsT, const double deltaT, const double cP, const double cS, const int indTemp, const double maxLen)
{
    int i, j;

    double vertsS[3][3];
    double pointF[3]; 
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            vertsS[i][j] = vertsT[9*blockIdx.x + 3*j + i];

    for (i = 0; i < 3; i++)
        pointF[i] = nodesMesh[3*blockIdx.y + i];
        

    double distMax = 0;
    double distMin = 100; //Mettere valore vicino a +Inf;
    double distCurr;

    double vettDist[3];

    for (i = 0; i < 3; i++)
    {
        vettDist[0] = vertsS[i][0] - pointF[0];
        vettDist[1] = vertsS[i][1] - pointF[1];
        vettDist[2] = vertsS[i][2] - pointF[2];
        distCurr = sqrt(vettDist[0]*vettDist[0] + vettDist[1]*vettDist[1] + vettDist[2]*vettDist[2]);
        if (distCurr > distMax)
            distMax = distCurr;
        if (distCurr < distMin)
            distMin = distCurr;
    }

    return ((indTemp - 2) * cS * deltaT > distMax + maxLen) || ((indTemp + 1) * cP * deltaT < distMin - 2*maxLen);
}


__device__ double dotProd3D(const double vettA[3], const double vettB[3])
{
    return vettA[0]*vettB[0] + vettA[1]*vettB[1] + vettA[2]*vettB[2];
}

__device__ double baseFunctionSM(const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM)
{
    double leftVector[3];
    leftVector[0] = vetCoeff[0] + matCoeff[0][0] * nodeInt[0] + matCoeff[0][1] * nodeInt[1] + matCoeff[0][2] * nodeInt[2];
    leftVector[1] = vetCoeff[1] + matCoeff[1][0] * nodeInt[0] + matCoeff[1][1] * nodeInt[1] + matCoeff[1][2] * nodeInt[2];
    leftVector[2] = vetCoeff[2] + matCoeff[2][0] * nodeInt[0] + matCoeff[2][1] * nodeInt[1] + matCoeff[2][2] * nodeInt[2];

    const double rightVector[3] = {indSM == 1, indSM == 2, indSM == 3};
    return dotProd3D(leftVector, rightVector);
}

__device__ void vectorRSM(double vectorRSMValue[3], const double normInt[3], const double matCoeff[3][3], const int indSM)
{
    vectorRSMValue[0] = + normInt[1] * matCoeff[indSM][2] - normInt[2] * matCoeff[indSM][1];
    vectorRSMValue[1] = - normInt[0] * matCoeff[indSM][2] + normInt[2] * matCoeff[indSM][0];
    vectorRSMValue[2] = + normInt[0] * matCoeff[indSM][1] - normInt[1] * matCoeff[indSM][0];
}


__device__ void nucleoKLT(double nuKLT[3][3], const double x[3], const double n[3], const double r, const double t, const double nodeInt[3], const double matCoeff[3][3], const double vetCoeff[3], const int indSM, 
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
            nuKLT[i][j] += nucleoTemp[i][j] * baseFunctionValue;
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


__device__ void nucleoKR(double nuKR[3][3], const double x[3], const double r, const double t, const double cP, const double cS, const double lambda, const double mu, const double rho, const int tensorE[3][3][3],
                                        const double normInt[3], const double matCoeff[3][3], const int indSM)
{
    int j, i, k;
    double nucleoCurr[3][3];
    double vettRSM[3];
    // Inizializzazione nucleo temporaneo
    for (i = 0; i < 3; i++)
        for (k = 0; k < 3; k++) 
            nucleoCurr[i][k] = 0;

    // Calcolo vettore Vsm
    vectorRSM(vettRSM, normInt, matCoeff, indSM - 1);

    // Ciclo sull'indice j dei nuclei nuKRj
    for(j = 0; j < 3; j++)
    {   
         // Calcolo nucleo corrente
         nucleoKRj(nucleoCurr, x, r, t, cP, cS, lambda, mu, rho, tensorE, j);
     
         // Ciclo sulle 9 componenti
         for (i = 0; i < 3; i++)
             for (k = 0; k < 3; k++) 
                 nuKR[i][k] -= vettRSM[j] * nucleoCurr[i][k];
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