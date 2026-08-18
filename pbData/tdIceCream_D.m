function g = tdIceCream_D(pbParam)
%TDICECREAM_D  Dirichlet datum loaded for the "tdIceCream" benchmark problem.
%   G = TDICECREAM_D(PBPARAM) returns a function handle G = @(x, t) -> (3x1 double) that
%   prescribes a time-oscillating displacement on the z-component only (zero on x/y),
%   with amplitude modulated by a 5th-power sine of t*(1 - (x1^2 + x2^2)/R^2) with 
%   R = 0.5, i.e. a radially-varying, ramping oscillation centered on the domain's axis.
%
%   Input arguments:
%       PBPARAM - (struct) physical/problem parameters (unused by this
%                 particular datum, but required by the GETDATUMHANDLEDIRICHLET
%                 calling convention), see READINPUTFILE.
%
%   Output arguments:
%       G - (function_handle) G(x, t) -> (3x1 double) displacement datum.
%
%   See also eebem.core.getDatumHandleDirichlet

arguments (Input)
    pbParam (1, 1) struct %#ok<INUSA>
end

arguments (Output)
    g (1, 1) function_handle
end

R = 0.5;
g = @(x, t) [0; 0; sin(t * (1 - ((x(1)^2 + x(2)^2) / R))).^5];

end