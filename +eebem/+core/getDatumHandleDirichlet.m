function g = getDatumHandleDirichlet(pbParam, basePath)
%GETDATUMHANDLEDIRICHLET  Return the exact Dirichlet (displacement) datum as a function handle.
%   G = GETDATUMHANDLEDIRICHLET(PBPARAM, BASEPATH) returns a function handle G = @(x, t)
%   with X a 1x3 point and T a scalar time, returning the 3x1 prescribed displacement.
%   For the built-in benchmark problems ("barH1", "barH3", "DesCop-cube",
%   "DesCop-sphere", "sphereWave", selected via PBPARAM.domainName) the closed-form
%   analytical solution is hard-coded. For any other PBPARAM.domainName, a user-supplied
%   "<domainName>_D.m" file is looked up under BASEPATH/pbData, loaded with
%   STR2FUNC/FEVAL, and expected to return a G = @(x, t) handle itself given PBPARAM.
%
%   Input arguments:
%       PBPARAM  - (struct) physical/problem parameters (domainName,
%                  velP, mu, nu, rho, Tfin, ...), see READINPUTFILE.
%       BASEPATH - (string, optional, default ".") project root, used
%                  to locate custom "pbData/<domainName>_D.m" files.
%
%   Output arguments:
%       G - (function_handle) G(x, t) -> (3x1 double) displacement datum.
%
%   Notes:
%       Asserts if PBPARAM.domainName is not one of the built-in
%       cases and no matching "pbData/<domainName>_D.m" file exists.
%
%   See also GETDATUMHANDLENEUMANN, CALCBOUNDDATADIRICHLET, CALCBETAI
arguments (Input)
    pbParam  (1, 1) struct
    basePath (1, 1) string = "."
end

arguments (Output)
    g function_handle
end

switch pbParam.domainName
    case {"barH1", "barH3"}
        h = 1 * strcmp(pbParam.domainName, "barH1") + 3 * strcmp(pbParam.domainName, "barH3");
        cP = pbParam.velP;
        hat_k = ceil((cP * pbParam.Tfin) / (2 .* h)) - 1;
        k_val = (0 : hat_k)';
        sol_an = @(x, t) sum((-1).^k_val .* ( (cP.*t - 2.*h.*k_val - (h-x(3))) .* ((cP.*t - 2.*h.*k_val - (h - x(3))) > 0) ...
                                        - (cP.*t - 2.*h.*(k_val+1) + (h-x(3))) .* ((cP.*t - 2.*h.*(k_val+1) + (h-x(3))) > 0)));
        g = @(x, t) [0; 0; sol_an(x, t)];

    case "DesCop-cube"
        a = 0.1;
        b = 100;
        cP = pbParam.velP;
        cS = pbParam.velS;
        rho = pbParam.rho;

        qP = @(r, t) sqrt(a) * (a*b - t + (r / cP));
        qS = @(r, t) sqrt(a) * (a*b - t + (r / cS));
        F = @(r, t) (1 / (2 * a * (r^2))) * (exp(-(qP(r, t)^2)) - exp(-(qS(r, t)^2)) + sqrt(a * pi) * (a*b - t) * (erf(qP(r, t)) - erf(qS(r, t))));
        
        fP = @(r, t) exp(-a * ((t - (r / cP) - (a*b))^2));
        fS = @(r, t) exp(-a * ((t - (r / cS) - (a*b))^2));

        d = @(i, j) (i == j);
        x0 = [1, 1, 1];

        u_ij = @(r, t, i, j) ((3 * r(i) * r(j) / (norm(r)^3)) - (d(i, j) / norm(r))) * F(norm(r), t) ...
                            + (r(i) * r(j) / (norm(r)^3)) * ((fP(norm(r), t) / (cP^2)) - (fS(norm(r), t) / (cS^2))) ...
                            + (d(i, j) / norm(r)) * (fS(norm(r), t) / (cS^2));

        u_i = @(r, t, i) u_ij(r, t, i, 1) + u_ij(r, t, i, 2) + u_ij(r, t, i, 3);
        
        u = @(r, t) (1 / (4*pi*rho)) .* [u_i(r, t, 1); u_i(r, t, 2); u_i(r, t, 3)];

        g = @(x, t) u(x - x0, t);

    case "DesCop-sphere"
        nu = pbParam.nu;
        cP = pbParam.velP;
        rho = pbParam.rho;
        R = 1;
        p0 = 143.5;
        a = 1279;
        b = 12792;

        omega0 = (cP * sqrt(1 - 2*nu)) / (R * (1 - nu));
        alpha0 = (cP * (1 - 2*nu)) / (R * (1 - nu));
        A2 = @(z) omega0 + (alpha0 - z)^2;
        A = @(z) sqrt(A2(z));
        beta0 = @(z) asin((alpha0 - z) / A(z));

       
        g_func = @(t, z) (p0 / (rho * A2(z) * omega0)) * ...
            ((1 / cP) * (-z * omega0 * exp(-z*t) + alpha0 * A(z) * exp(-alpha0*t) * cos(omega0*t - beta0(z)) + omega0 *  A(z) * exp(-alpha0*t) * sin(omega0*t - beta0(z))) ...
              - (1 / R) * (A(z) * exp(-alpha0*t) * cos(omega0*t - beta0(z)) - omega0 * exp(-z*t)));

        g = @(x, t) (g_func(t, a) - g_func(t, b)) .* (-x' ./ norm(x)); 

    case "sphereWave"
        cP = pbParam.velP;
        ondaIncid = @(x, t) exp(-20 .* ((x(1) - 2 + cP.*t - 0.475).^2));
        g = @(x, t) [ondaIncid(x, t); 0; 0];

    otherwise
        fileName = pbParam.domainName + "_D.m";
        assert(exist(fullfile(basePath, "pbData", fileName), 'file'), "Datum file not found. Provide a .m file returning the necessary function handle");
        func = str2func(extractBefore(fileName, "."));
        addpath(fullfile(basePath, "pbData"))
        g = feval(func, pbParam);
        rmpath(fullfile(basePath, "pbData"))
end
end