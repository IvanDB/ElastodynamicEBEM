function g = industrialComponent_N(pbParam)
%INDUSTRIALCOMPONENT_N  Neumann datum for the "industrialComponent" benchmark problem.
%   G = INDUSTRIALCOMPONENT_N(PBPARAM) returns a function handle G = @(x, t, n) -> (3x1
%   double) that prescribes a constant-in-time, piecewise-constant-in-space traction:
%   zero on the first two components, and +/-1 on the third (sign matching the outward
%   normal N's z-component) wherever |N(3)| > 0.5, i.e. on roughly horizontal facets of
%   the mesh, and 0 elsewhere. Loaded dynamically by GETDATUMHANDLENEUMANN via
%   "pbData/industrialComponent_N.m" when PBPARAM.domainName == "industrialComponent".
%
%   Input arguments:
%       PBPARAM - (struct) physical/problem parameters (unused by this
%                 particular datum, but required by the GETDATUMHANDLENEUMANN
%                 calling convention), see READINPUTFILE.
%
%   Output arguments:
%       G - (function_handle) G(x, t, n) -> (3x1 double) traction datum.
%
%   See also GETDATUMHANDLENEUMANN

arguments
    pbParam (1, 1) struct
end

g = @(x, t, n) [0; 0; sign(n(3)) * (abs(n(3)) > 0.5)];
end
