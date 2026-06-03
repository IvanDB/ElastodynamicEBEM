function matrixW = calcMatrixW(matrixSpecs, nGPU, basePath, pbParam, domainMesh, quadData, constValues)
%CALCMATRIXV Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    matrixSpecs struct
    nGPU        (1, 1) double {mustBeInteger, mustBePositive}
    basePath    (1, 1) string
    pbParam     struct
    domainMesh  struct
    quadData    struct
    constValues cell
end

arguments (Output)
    matrixW cell
end

import eebem.core.*
import eebem.utility.*

%Allocazione array contente i blocchi matriciali
matrixW = cell(matrixSpecs.maxNumBlocksPerIter, matrixSpecs.numIter);

for indIter = 1 : matrixSpecs.numIter
    %Calcolo numero blocchi iterazione corrente
    numBlocksThisIter = matrixSpecs.offsets_full(indIter + 1) - matrixSpecs.offsets_full(indIter);

    %Avvio computazione GPU
    % spmd(nGPU)
        indGPU = 1; %spmdIndex;

        globIdx = nGPU * (indIter - 1) + indGPU;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMulti = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(indGPU); % gpuDevice(spmdIndex);

            srcPath = fullfile(basePath, "+eebem", "+core", "kernelsCUDA", "kernelW.cu");
            ptxPath = fullfile(basePath, "buildDir", "kernelW.ptx");

            kernelW = parallel.gpu.CUDAKernel(ptxPath, srcPath);
            kernelW.GridSize = [];
            kernelW.ThreadBlockSize = [quadData.methodSpecs.numSRint quadData.methodSpecs.numGHint 1];
            kernelW.SharedMemorySize = quadData.methodSpecs.numINT * 9 * 8;

            gpuInputArrays = copyArrayW(domainMesh, quadData, constValues);

            matrixOutMulti = launchCUDAKernel(gpuID, kernelW, pbParam.deltaT, pbParam.velP, pbParam.velS, pbParam.lambda, pbParam.mu, pbParam.rho, 4*pi*pbParam.deltaT, ...
                                                gpuInputArrays.stdEXTw, gpuInputArrays.stdEXTnx, gpuInputArrays.stdEXTny, gpuInputArrays.stdEXTnz, quadData.methodSpecs.numEXT, ...
                                                gpuInputArrays.stdINTw, gpuInputArrays.stdINTnx, gpuInputArrays.stdINTny, gpuInputArrays.stdINTnz, ...
                                                gpuInputArrays.vertsT, gpuInputArrays.areeT, gpuInputArrays.normT, gpuInputArrays.indSMmatrix, gpuInputArrays.matCoeff, gpuInputArrays.vetCoeff, ...
                                                matrixSpecs.offsets_sing(globIdx), matrixSpecs.numBlocks, gpuInputArrays.nodesMesh, domainMesh.lMax);
        end
    %end

    %Allocazione array contente i componenti singolari dei blocchi di questa iterazione
    % matrixSubBlocksSING = zeros(3,3,?, matrixSpecs.blockSizes2D(1), numBlocksThisIter);
    matrixSubBlocksSING = cell(matrixSpecs.blockSizes2D, numBlocksThisIter);% Cell array of same dimentrions as the 
    % GPU output. Each cell will represent the interaction between two
    % nodes at a specific time instant (the same aa a GPU block)and so it
    % wil be a 3x3 matrix. To compute it we use the same function we
    % already build but instead of skipping ovelapping triangles we skip
    % the ones that are different. In that way that 3x3 block represents
    % the exact contribution of the singular itegrations to be added to the
    % exact same block of the GPU matrix. (we just need to sum the two matrixes after cell2mat)

    % The first two dimentions represent the actual matrix block. The third
    % dimention is the maximum amount of singular contributions that an
    % outer node can give (itself plus al the triangles that touch it)

    %Avvio computazione CPU
    % Need to structure this part and final matrix assembly
end
end









function gpuInputArrays = copyArrayW(domainMesh, quadData, constValues)

end