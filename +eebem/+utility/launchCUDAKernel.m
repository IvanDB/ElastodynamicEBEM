function deviceMatrix = launchCUDAKernel(gpuID, kernel, varargin)
%LAUNCHCUDAKERNEL  Allocate the output buffer and run a compiled CUDA kernel on one GPU.
%   DEVICEMATRIX = LAUNCHCUDAKERNEL(GPUID, KERNEL, VARARGIN) waits for GPUID to be
%   idle, allocates a zero gpuArray output buffer sized from KERNEL.GridSize (9
%   vector components per grid cell, matching the 3x3 block layout used throughout
%   this codebase), and invokes KERNEL with that buffer as first argument followed
%   by VARARGIN. Used by CALCMATRIXV, CALCMATRIXK, CALCMATRIXK_C and CALCMATRIXW
%   to launch their respective "kernel*.cu" CUDA kernels.
%
%   Input arguments:
%       GPUID    - (parallel.gpu.GPUDevice) the target GPU device.
%       KERNEL   - (parallel.gpu.CUDAKernel) the compiled kernel to
%                  run, with GridSize/ThreadBlockSize/SharedMemorySize
%                  already configured by the caller.
%       VARARGIN - (repeating) extra arguments forwarded to KERNEL after
%                  the output buffer (typically physical constants
%                  followed by quadrature/mesh gpuArray inputs).
%
%   Output arguments:
%       DEVICEMATRIX - (gpuArray double) the kernel's output
%                      buffer, still on the GPU (not gathered).
%
%   See also CALCMATRIXV, CALCMATRIXK, CALCMATRIXK_C, CALCMATRIXW
arguments
    gpuID (1,1) parallel.gpu.GPUDevice
    kernel (1,1) parallel.gpu.CUDAKernel
end

arguments (Repeating)
    varargin
end

%Array di input
linearMatrixSixe = 9 * prod(kernel.GridSize);

deviceMatrix = zeros(linearMatrixSixe, 1, 'gpuArray');
wait(gpuID);

%Avvio computazione GPU
deviceMatrix = feval(kernel, deviceMatrix, varargin{:});

end