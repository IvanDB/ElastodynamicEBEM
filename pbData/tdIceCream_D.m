function g = ConeBall_D(pbParam)
%CONEBALL_D  Dirichlet datum loaded for the "tdIceCream" benchmark problem.
%   G = CONEBALL_D(PBPARAM) returns a function handle G = @(x, t) -> (3x1 double) that
%   prescribes a time-oscillating displacement on the z-component only (zero on x/y),
%   with amplitude modulated by a 5th-power sine of t*(1 - (x1^2+x2^2)/R^2) and R =
%   0.5, i.e. a radially-varying, ramping oscillation centered on the domain's axis.
%
%   Input arguments:
%       PBPARAM - (struct) physical/problem parameters (unused by this
%                 particular datum, but required by the GETDATUMHANDLEDIRICHLET
%                 calling convention), see READINPUTFILE.
%
%   Output arguments:
%       G - (function_handle) G(x, t) -> (3x1 double) displacement datum.
%
%   Notes:
%       NAME MISMATCH: this file is "tdIceCream_D.m" (and is loaded by
%       GETDATUMHANDLEDIRICHLET via STR2FUNC("tdIceCream_D") when PBPARAM.domainName ==
%       "tdIceCream"), but the function defined inside it is named "ConeBall_D", not
%       "tdIceCream_D". MATLAB resolves the callable name from the file name, so
%       STR2FUNC/FEVAL("tdIceCream_D", ...) does work here -- but the function's own name
%       (as seen by e.g. NARGIN, error stacks, or a direct "help ConeBall_D") does not
%       match the file, which is a likely leftover from copying/renaming this datum file
%       from an (unseen) "ConeBall_D.m" original. Worth reconciling one way or the other.
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