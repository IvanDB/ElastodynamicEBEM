function OutputMatrix = calcSingSubBlockW(pbParam, domainMesh, quadData, constData, indSMatrix, TriangPerNodes, maxNumTriangles, timeInstant, outerNode, innerNode)
%CALCSINGSUBBLOCKW  Correct the singular vertex-pair contribution of a hypersingular (W) matrix block.
%   OUTPUTMATRIX = CALCSINGSUBBLOCKW(PBPARAM, DOMAINMESH, QUADDATA, CONSTDATA, INDSMATRIX,
%                                    TRIANGPERNODES, MAXNUMTRIANGLES, TIMEINSTANT, OUTERNODE, INNERNODE)
%   evaluates the singular correction to add between mesh vertices OUTERNODE and
%   INNERNODE (only meaningful, and only called by CALCMATRIXW, for vertex pairs
%   that share at least one triangle). For every shared triangle it delegates to
%   the local helper COMPUTESINGBLOCKW, which evaluates the closed-form regularized
%   hypersingular kernel (see the local KERNELCALC family) on the
%   light-cone-intersected sub-triangles and combines it with a 4-point backward
%   finite difference in time (coefficients [1, -3, 3, -1], scaled by 1/deltaT^2)
%   to realize the double time-differentiation the hypersingular operator requires.
%
%   Input arguments:
%       PBPARAM         - (struct) physical/time-discretization
%                         parameters, see READINPUTFILE.
%       DOMAINMESH      - (struct) triangulated boundary mesh (normal), see READSPACEMESH.
%       QUADDATA        - (struct) quadrature nodes/weights/METHODSPECS,
%                         see GENERATEQUADDATA.
%       CONSTDATA       - (cell) per-triangle data from CALCCONSTVALUES.
%       INDSMATRIX      - (numVertices x maxNumTriangles int32) for each vertex
%                         and each incident-triangle slot, the local index (1, 2
%                         or 3) of that vertex within the triangle (0 if the slot
%                         is unused padding); built locally by the COPYARRAYW
%                         helper in CALCMATRIXW, not by READSPACEMESH.
%       TRIANGPERNODES  - (numVertices x maxNumTriangles) triangles
%                         incident to each vertex (0-padded).
%       MAXNUMTRIANGLES - (positive integer) width of TRIANGPERNODES,
%                         i.e. the maximum vertex valence in the mesh.
%       TIMEINSTANT     - (nonnegative integer) discrete time-lag index for this block.
%       OUTERNODE       - (positive integer) global index of the test vertex.
%       INNERNODE       - (positive integer) global index of the trial vertex.
%
%   Output arguments:
%       OUTPUTMATRIX - (3x3 double) singular correction to add at the
%                      (OUTERNODE, INNERNODE) position of the block.
%
%   See also CALCMATRIXW, CALCCONSTVALUES, GENERATEFINALG2DNODES

arguments (Input)
    pbParam         (1, 1) struct
    domainMesh      (1, 1) struct
    quadData        (1, 1) struct
    constData       (:, 1) cell
    indSMatrix      (:, :) int32
    TriangPerNodes  (:, :) int32
    maxNumTriangles (1, 1) double {mustBeInteger, mustBePositive}
    timeInstant     (1, 1) double {mustBeInteger, mustBeNonnegative}
    outerNode       (1, 1) double {mustBeInteger, mustBePositive}
    innerNode       (1, 1) double {mustBeInteger, mustBePositive}
end

arguments (Output)
    OutputMatrix (3, 3) double
end
delta = eye(3);

epsilon = zeros(3, 3, 3); 


epsilon(1, 2, 3) = 1;
epsilon(2, 3, 1) = 1;       
epsilon(3, 1, 2) = 1;


epsilon(1, 3, 2) = -1;
epsilon(2, 1, 3) = -1;
epsilon(3, 2, 1) = -1;


OutputMatrix = zeros(3,3);

outerTriangles = TriangPerNodes(outerNode, :); 
innerTriangles = TriangPerNodes(innerNode, :);
outerVertexes = indSMatrix(outerNode,:); 
                                                
innerVertexes = indSMatrix(innerNode,:);


for i = 1 : maxNumTriangles 
    currentOuterTriangle = outerTriangles(i); 
    if (currentOuterTriangle == 0)
        break;
    end
    currentOuterVertex = outerVertexes(i);
    for j = 1 : maxNumTriangles 
        currentInnerTriangle = innerTriangles(j); 
        if (currentInnerTriangle == 0)
        break;
        end
        if (currentInnerTriangle ~= currentOuterTriangle)
            % if the two triangles aren't equal it does nothing and then
            % proceeds
            continue;
        end

        currentInnerVertex = innerVertexes(j);
        
        currentOuterNormal = domainMesh.normal(currentOuterTriangle,:); 

        singularSubBlock = computeSingBlockW(pbParam, currentOuterVertex, currentInnerVertex, currentOuterNormal, constData{currentOuterTriangle}, timeInstant, quadData.methodSpecs, delta, epsilon);
        OutputMatrix = OutputMatrix + singularSubBlock; 
    end
end

end



function VXi = VFunction(triangleCoeffs, triangleMatrix, vertexInd, xi)
%Evaluate the piecewise-linear trial basis function of local vertex
%VERTEXIND (the "V_alpha^s(xi)" of the regularized kernel) at point XI.
eVec = [0,0,0];
eVec(vertexInd) = 1;
VXi = dot((triangleCoeffs + triangleMatrix*xi'), eVec); %traspongo xi perché viene passato come riga
end

function VTildeX = VTildeFunction(triangleCoeffs, triangleMatrix, vertexInd, x)
%Evaluate the piecewise-linear test basis function of local vertex VERTEXIND
%(the "V_tilde-alpha^tilde-s(x)" of the regularized kernel) at point X.
eVec = [0,0,0];
eVec(vertexInd) = 1;
VTildeX = dot((triangleCoeffs +triangleMatrix*x'), eVec); %traspongo x perché viene passato come riga
end

function kernelTot = kernelCalc(t, VTildeX, VXi, VAlphaS, VTildeAlphaS, rVector, r, n, v, mu, lambda, rho, cS, cP, delta, epsilon)
%Regularized hypersingular kernel at one integration point: sum of the
%"T" (tangential/basis-weighted) and "R" (radial/curl-weighted) parts.
kernelTot = zeros(3,3);
kernelT = zeros(3,3);
kernelR = zeros(3,3);
for i = 1 : 3
    for k = 1 : 3
        kernelT(i,k) = kernelTCalc(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, rho, cS, cP, delta);
        kernelR(i,k) = kernelRCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon);
        kernelTot(i,k)= kernelR(i,k) + kernelT(i,k);
    end
end
end


function kernelT = kernelTCalc(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, rho, cS, cP, delta)
%"T" part of the kernel: sum of two Heaviside-jump ("Delta1"/"Delta2")
%terms and one smooth Heaviside-integrated ("H") term.

kernelTDelta1 = calcKernelTDelta1(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta);
kernelTDelta2 = calcKernelTDelta2(i, k, t, VTildeX, VXi, r, n, v, mu, lambda, rho, cS, cP, delta);
kernelTH = calcKernelTH(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta);

kernelT = kernelTDelta1 + kernelTDelta2 + kernelTH;
end


function kernelTDelta1 = calcKernelTDelta1(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta)
%First (1/r) Heaviside-jump term of the "T" kernel part.
H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IDeltaS1 = H_S/(cS*cS);
IDeltaP1 = H_P/(cP*cP);

Rv = dot(rVector,v);
Rn = dot(rVector,n); % Notazione per prodotti scalari
Nv = dot(n,v);

constTerm = 1/(2*pi*(lambda+mu));

kernelTDelta1 = VTildeX*VXi*constTerm*(lambda*lambda*((n(k)*v(i))/(2*r)) + lambda*mu*((rVector(i)*n(k)*Rv + rVector(k)*v(i)*Rn)/(r*r*r)) +...
                    mu*mu*((rVector(k)*n(i)*Rv + rVector(i)*v(k)*Rn + delta(i,k)*Rv*Rn + rVector(i)*rVector(k)*Nv)/(2*r*r*r)))*(IDeltaS1-IDeltaP1);
end

function kernelTDelta2 = calcKernelTDelta2(i, k, t, VTildeX, VXi, r, n, v, mu, lambda, rho, cS, cP, delta)
%Second (1/r, higher wave-speed power) Heaviside-jump term of the "T" kernel part.
H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IDeltaS2 = H_S/(cS*cS*cS*cS); 
IDeltaP2 = H_P/(cP*cP*cP*cP);

Nv = dot(n,v);

constTerm = -(lambda*mu + 2*mu*mu)/(4*pi*rho*(lambda+mu));

kernelTDelta2 = VTildeX*VXi*constTerm*((lambda*n(k)*v(i) + mu*n(i)*v(k) + mu*delta(i,k)*(Nv))/(r))*(IDeltaS2-IDeltaP2);
end

function kernelTH = calcKernelTH(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta)
%Smooth, Heaviside-integrated ("H", i.e. time-integrated twice) term of the "T" kernel part.

H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IHS = H_S*(((t^2)-(r/cS)^2)/2);
IHP = H_P*(((t^2)-(r/cP)^2)/2);

Rn = dot(rVector,n);
Rv = dot(rVector,v);
Nv = dot(n,v);

constTerm = mu/(2*pi*(lambda+mu));

kernelTH = VTildeX*VXi*constTerm*(3*lambda*((rVector(k)*v(i)*Rn + rVector(i)*n(k)*Rv)/(r*r*r*r*r)) - 2*lambda*((n(k)*v(i))/(r*r*r))...
               +3*mu*((rVector(k)*n(i)*Rv + rVector(i)*v(k)*Rn + delta(i,k)*Rv*Rn + rVector(i)*rVector(k)*Nv)/(2*r*r*r*r*r)) - mu*((n(i)*v(k) + delta(i,k)*Nv)/(r*r*r)))*(IHS-IHP);
end


function kernelR = kernelRCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon)
%"R" part of the kernel: sum of a Heaviside-jump ("Delta") and a smooth Heaviside-integrated
%("H") term, both built from the curl-weighted vectors VALPHAS/VTILDEALPHAS.

kernelRDelta = calcKernelRDelta(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta);
kernelRH = calcKernelRH(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta);

kernelR = kernelRH + kernelRDelta;

end

function kernelRH = calcKernelRH(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta)
    %Smooth, Heaviside-integrated ("H") term of the "R" kernel part.
    H_P = (t - (r/cP) > 0);
    H_S = (t - (r/cS) > 0);
    IHRGS = H_S*(((t-(r/cS))^4)/(24) + r*((t-(r/cS))^3)/(6*cS));
    IHRGP = H_P*(((t-(r/cP))^4)/(24) + r*((t-(r/cP))^3)/(6*cP));
    
    % --- PRE-CALCOLO VETTORIALE ---
    A = cross(rVector, VTildeAlphaS);
    B = cross(rVector, VAlphaS);
    Vdot = dot(VTildeAlphaS, VAlphaS);
    
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    
    % --- FORMULA CONTRATTA ---
    % Termini divisi per r^5
    term_r5 = (3/r5) * ( 2*lambda * A(i) * B(k) ...
                       + (lambda + 2*mu) * B(i) * A(k) ...
                       - (lambda + 2*mu) * Vdot * rVector(i) * rVector(k) );
                   
    % Termini divisi per r^3
    term_r3 = (1/r3) * ( -2*lambda * delta(i,k) * Vdot ...
                       + 2*lambda * VAlphaS(i) * VTildeAlphaS(k) ...
                       + (lambda + 2*mu) * VTildeAlphaS(i) * VAlphaS(k) );
                   
    valoreContratto = term_r5 + term_r3;
    
    constTerm = -((mu*mu)/(4*pi*rho*(lambda+mu)));
    kernelRH = constTerm * valoreContratto * (IHRGS - IHRGP);
end

function kernelRDelta = calcKernelRDelta(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta)
    %Heaviside-jump ("Delta") term of the "R" kernel part.
    H_P = (t - (r/cP) > 0); 
    H_S = (t - (r/cS) > 0);
    IDeltaRGS = H_S*(((t-(r/cS))^2)/(2*cS*cS));
    IDeltaRGP = H_P*(((t-(r/cP))^2)/(2*cP*cP));
    
    % --- PRE-CALCOLO VETTORIALE (Sostituisce i 6 cicli!) ---
    A = cross(rVector, VTildeAlphaS); % Prodotto vettoriale
    B = cross(rVector, VAlphaS);      % Prodotto vettoriale
    Vdot = dot(VTildeAlphaS, VAlphaS); % Prodotto scalare (V_Tilde * V)
    
    r2 = r*r;
    r3 = r2*r;
    
    % --- FORMULA CONTRATTA ---
    valoreContratto = (2*lambda / r3) * A(i) * B(k) ...
                    + ((lambda + 2*mu) / r3) * B(i) * A(k) ...
                    + (lambda + 2*mu) * Vdot * (delta(i,k)/r - (rVector(i)*rVector(k))/r3);
    
    constTerm = -((mu*mu)/(4*pi*rho*(lambda+mu)));
    kernelRDelta = constTerm * valoreContratto * (IDeltaRGS - IDeltaRGP);
end

function singularSubBlock = computeSingBlockW(pbParam, outerVertex, innerVertex, triangleNormal, currentConstData, timeInstant, methodSpecs, delta, epsilon)
%Integrate the regularized hypersingular kernel over one shared triangle (area integral via
%light-cone sub-triangles) and combine it with the 4-point backward time difference, for
%one (outerVertex, innerVertex) pair. Called once per shared triangle by CALCSINGSUBBLOCKW.

singularSubBlock = zeros(3,3);

mu = pbParam.mu;
lambda = pbParam.lambda;
rho = pbParam.rho;
deltaT = pbParam.deltaT;
cS = pbParam.velS;
cP = pbParam.velP;

kernelCoeffs = [1, -3, 3, -1]; % I coefficienti C_\eta
timeCoeffs = [-2, -1, 0, 1]; % Gli \eta che individuano i vari istanti temporali
temporalInstants = timeInstant + timeCoeffs;

numExtSubRegion = methodSpecs.numSRext;
numExtNodesPerSubRegion = methodSpecs.numGHext;

minDiffTemp = temporalInstants(1)*deltaT;

triangleCoeffs = currentConstData.vetCoeff;
triangleMatrix = currentConstData.matCoeff;

OuterVVector = cross(triangleNormal,triangleMatrix(outerVertex,:));
InnerVVector = cross(triangleNormal,triangleMatrix(innerVertex,:));

for k = 1 : 4 % Ciclo su istanti temporali
    if (temporalInstants(k)<=0)
        continue;
    end
    currentTime = temporalInstants(k)*deltaT;
    currentCoeff = kernelCoeffs(k);
    outerIntegral = zeros(3,3);

    for indSR = 1 : numExtSubRegion %ciclo su sottoregioni e nodi
        for indNode = 1 : numExtNodesPerSubRegion

            indEXTn = (indSR-1) * methodSpecs.numGHext + indNode; %Trovo l'indice globale del nodo esterno corrente
            currentOuterNode = currentConstData.EXTn{indEXTn};
            currentOuterWeight = currentConstData.EXTw{indNode};

            VTildeX = VTildeFunction(triangleCoeffs, triangleMatrix, outerVertex, currentOuterNode);

            innerIntegral = zeros(3, 3);

            for indChild = 1 : 3 %ciclo sui triangoli figli interni
                rMin = max(cS * minDiffTemp, 0);
                rInt = cS * currentTime;
                rExt = cP * currentTime;

                [G2DCnodes, G2DCweights] = eebem.utility.quadratureRules.generateFinalG2Dnodes(currentConstData.childVerts{indEXTn, indChild}, rMin, rInt, rExt, numNodes = sqrt(methodSpecs.numSNGLR));
                
                for indInnerNode = 1 : length(G2DCnodes) % ciclo sui nodi interno del figlio corrente
                    currentInnerNode = G2DCnodes(indInnerNode,:);
                    currentInnerWeight = G2DCweights(indInnerNode);

                    VXi = VFunction(triangleCoeffs, triangleMatrix, innerVertex, currentInnerNode);

                    rVector = currentOuterNode - currentInnerNode;
                    rNorm = norm(rVector,2);

                    kernel = kernelCalc(currentTime, VTildeX, VXi, InnerVVector,OuterVVector,...
                        rVector, rNorm, triangleNormal, triangleNormal, mu, lambda, rho, cS, cP, delta, epsilon);

                    innerIntegral = innerIntegral + currentInnerWeight*kernel;
                end
            end
            outerIntegral = outerIntegral + currentOuterWeight*innerIntegral;
        end
    end
    singularSubBlock = singularSubBlock + (currentCoeff/(deltaT^2))*outerIntegral;
end
end




