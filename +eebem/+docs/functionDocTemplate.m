function output = templateFunctionName(input1, input2, name1, name2)
%TEMPLATEFUNCTIONNAME  One-line summary, imperative mood, under ~75 chars.
%   OUTPUT = TEMPLATEFUNCTIONNAME(INPUT1, INPUT2) explains in one or two
%   sentences what the function computes or does. Refer to input/output
%   names in UPPERCASE, exactly as they appear in the function signature,
%   so the text stays consistent if the code is skimmed alongside it.
%
%   OUTPUT = TEMPLATEFUNCTIONNAME(..., NAME1, VALUE1) documents an additional
%   calling form (e.g. a name-value pair). Add one such paragraph per meaningful
%   syntax; omit this block if the function has a single calling form.
%
%   Input arguments:
%       INPUT1 - (class, size) role of the argument in this function.
%       INPUT2 - (class, size) role of the argument in this function.
%       NAME1  - (class, size) optional name-value argument.
%                State the default value and what it controls.
%
%   Output arguments:
%       OUTPUT - (class, size) meaning of the returned value.
%
%   Notes:
%       Optional paragraph, keep only if genuinely useful. Typical uses in this project:
%       GPU/parallel-pool requirements, files read or written on disk, numerical/physical
%       assumptions, known limitations, "deprecated" or "work in progress" notices.
%
%   See also RELATEDFUNCTIONONE, RELATEDFUNCTIONTWO

arguments
    input1
    input2
    name1.someName = "defaultValue"
    name2.otherName (1, 1) double = 0
end

output = []; %#ok<NASGU> placeholder body, template only
end

%{
=========================================================================
 WHY THIS LAYOUT
=========================================================================
MATLAB's HELP/DOC/LOOKFOR only ever look at the block of contiguous "%"
lines that comes IMMEDIATELY after the "function ..." line. If a comment
appears after an "arguments" block instead, HELP does not show it. That
is why the doc block above sits BEFORE "arguments", matching the pattern
MathWorks' own style guide recommends and matching most of the
auto-generated stubs already present in this codebase.

  - H1 line (first line): starts with %FUNCTIONNAME in caps, no blank
    line before it. This is the line shown by `lookfor` and in folder
    listings, so it must be a self-contained one-line summary.
  - Blank "%" separators are fine (they do not end the help block, only
    a truly empty line / first line of code does).
  - "See also" is a plain comment line; MATLAB turns the names into
    clickable links automatically if those functions are on the path.

=========================================================================
 PROJECT-WIDE GLOSSARY (recurring struct/cell arguments)
=========================================================================
Several structs/cells are threaded through most of +eebem.+core and
+eebem.+utility unchanged. To keep individual doc blocks short, describe
them tersely and consistently using this glossary instead of re-deriving
their fields every time:

  basePath      (1,1) string   - Root folder of the repository; used to
                                  build paths to inputFiles/, mesh/,
                                  pbData/, buildDir/, tempData/,
                                  outputData/, outputPlot/.
  pbParam       (1,1) struct   - Physical + discretization parameters of
                                  the problem: rho, mu, lambda, nu, velP,
                                  velS, domainName, meshName, Tfin, nT,
                                  deltaT, beta, tMlt, STcoupling, BIE,
                                  BOU, defaultValues, ... (built by
                                  utility.fileRead.readInputFile and
                                  completed by utility.finalizeParameters).
  domainMesh    (1,1) struct   - Triangular boundary mesh of the scatterer
                                  domain: coordinates, triangles, normal,
                                  area, center, curl, maxL, indSMmatrix,
                                  numVertices, numTriangles, name, lev
                                  (built by utility.fileRead.readSpaceMesh).
  quadData      (1,1) struct   - Quadrature nodes/weights (external,
                                  internal, 1D singular) plus
                                  quadData.methodSpecs (built by
                                  utility.generateQuadData).
  methodSpecs   (1,1) struct   - Quadrature method configuration nested
                                  inside quadData (node counts, quadType,
                                  stringID).
  matrixSpecs   (1,1) struct   - GPU memory/iteration layout for a
                                  block-Toeplitz BEM matrix (block sizes,
                                  number of blocks, per-device offsets),
                                  built by core.calcMatrixSpecs.
  constValues   (:,1) cell     - Per-triangle precomputed geometric data
                                  (P1 basis coefficients, quadrature nodes
                                  mapped to the triangle, singular
                                  sub-triangle data), built by
                                  core.calcConstValues.
  fullFileNames (1,1) struct   - Output .mat file paths (tmpFullFilename,
                                  outFullFilename), built by
                                  utility.generateFilenames.

=========================================================================
 SCOPE OF THIS DOCUMENTATION PASS
=========================================================================
A full HELP-style block (as above) is added to the PRIMARY function of
every .m file (the one whose name matches the file name, which is the
only one HELP/LOOKFOR can ever address). Local/private helper functions
defined further down in the same file (e.g. copyArrayV in calcMatrixV.m,
or the kernel* helpers in calcSingSubBlockW.m) get a single-line "%"
comment instead: MATLAB never exposes their help text independently, so
a full block there would just add noise.
%}
