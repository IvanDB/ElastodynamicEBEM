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
    spmd(nGPU)
        indGPU = spmdIndex;
        globIdx = nGPU * (indIter - 1) + indGPU;
        numBlockThisLaunch = matrixSpecs.offsets_sing(globIdx + 1) - matrixSpecs.offsets_sing(globIdx);

        matrixOutMulti = [];
        
        if(numBlockThisLaunch > 0)
            gpuID = gpuDevice(spmdIndex);

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
    end

    %Allocazione array contente i componenti singolari dei blocchi di questa iterazione
    matrixSubBlocksSING = zeros(); % Need to find a way to organize this for kernelW

    %Avvio computazione CPU
    % Need to structure this part and final matrix assembly
end
end









function gpuInputArrays = copyArrayW(domainMesh, quadData, constValues)

end