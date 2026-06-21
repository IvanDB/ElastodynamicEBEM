function autobuild(basePath, flag)

arguments
    basePath (1, 1) string
    flag     (1, 1) logical = false
end

disp("Autobuilding is in beta")
addpath(fullfile(basePath, "buildDir"));

if(~flag)
    return
end

compilerInfo = mex.getCompilerConfigurations('C');
if(isempty(compilerInfo))
    mex -setup C
    compilerInfo = mex.getCompilerConfigurations('C');
end
assert(~isempty(compilerInfo), "C compiler can't be found")

scriptCompiler = "";
if(ispc)
    assert(contains(compilerInfo.ShortName, "MSVC"), "For CUDA compilation on Windows the host compiler must be MSVC")

    scriptCompiler = ['"' compilerInfo.Details.CommandLineShell '" ' ...
        compilerInfo.Details.CommandLineShellArg ' && '];
end


coreMEXCkernelsDirectory = fullfile(basePath, "+eebem", "+core", "kernelsMEXC");
coreMEXCkernelsFileNames = ["kernelK.c", "kernelV.c"];

binOutputDirectory = fullfile(basePath, "buildDir");

for MEXCkernel = fullfile(coreMEXCkernelsDirectory, coreMEXCkernelsFileNames)
    cmd_base = "mex -v -R2018a ";
    if(contains(compilerInfo.ShortName, "MSVC"))
        cmd_opts = "LINKFLAGS='$LINKFLAGS /LTCG' COMPFLAGS='$COMPFLAGS /O2 /GL /fp:fast' ";
    end
    if(contains(compilerInfo.ShortName, "gcc") || contains(compilerInfo.ShortName, "clang") || contains(compilerInfo.ShortName, "mingw"))
        cmd_opts = "LDFLAGS='$LDFLAGS -flto -march=native -mtune=native' CFLAGS='$CFLAGS -std=c99 -O3 -flto -ffast-math -march=native -mtune=native' ";
    end

    cmd_outf = "-outdir '" + binOutputDirectory + "' "; 

    cmd = cmd_base + cmd_opts + cmd_outf + " '" + MEXCkernel + "' ";
    cmdout = evalc(cmd);
    assert(contains(cmdout, "MEX completed successfully."), "C-MEX compilation failed! Shell output:" + newline + cmdout)
end


gpuInfo = gpuDevice;
gpuCC = gpuInfo.ComputeCapability;
gpuCC(gpuCC == '.') = [];

coreCUDAkernelsDirectory = fullfile(basePath, "+eebem", "+core", "kernelsCUDA");
coreCUDAkernelsFileNames = ["kernelK.cu", "kernelV.cu", "kernelKboundary.cu", "kernelKinternal.cu"];

for CUDAkernel = fullfile(coreCUDAkernelsDirectory, coreCUDAkernelsFileNames)
    cmd = strcat(scriptCompiler, ' nvcc -ptx -O3 -Wno-deprecated-gpu-targets -arch=compute_', gpuCC, ' -odir "', binOutputDirectory, '" "', CUDAkernel, '"');
    [status, cmdout] = system(cmd);
    assert(~status, "CUDA compilation failed! Shell output:" + newline + cmdout)
end
end