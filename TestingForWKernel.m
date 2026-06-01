function OutputMatrix = TestingForWKernel(outerNode, innerNode, timeInstant, indSMatrix,...
    TriangPerNodes, pbParam, domainMesh, quadData, constData, maxNumTriangles)


% test per calcolare un blocchetto di matrice del W
% si passano i nodi interni ed esterni (s e sTilde) come inner e outer node
% si passa l'istante temporale l come timeInstant
% questa funzione si occupa di calcolare il blocchetto 3x3
% [W_{\tlde{s},s}^{(l)}]_{ik} o almeno parti di essa. Così controlliamo se
% il nucleo è giusto

% arguments (Input)
%     pbParam
%     domainMesh
%     quadData
% end

arguments (Output)
    OutputMatrix
end
delta = eye(3); % delta di kronecher

epsilon = zeros(3, 3, 3); 

% Permutazioni pari (+1)
epsilon(1, 2, 3) = 1;
epsilon(2, 3, 1) = 1;       % Per simboli levi-Civita
epsilon(3, 1, 2) = 1;

% Permutazioni dispari (-1)
epsilon(1, 3, 2) = -1;
epsilon(2, 1, 3) = -1;
epsilon(3, 2, 1) = -1;


OutputMatrix = zeros(3,3);
mu = pbParam.mu;
lambda = pbParam.lambda;
rho = pbParam.rho;
deltaT = pbParam.Tfin/pbParam.nT;
cS = pbParam.velS;
cP = pbParam.velP;
kernelCoeffs = [1,-3,3,1]; % I coefficienti C_\eta
timeCoeffs = [-2,-1,0,1]; % Gli \eta che individuano i vari istanti temporali
outerTriangles = TriangPerNodes(outerNode, :); % Trovo gli indici dei triangoli che hanno outerNode come vertice 
innerTriangles = TriangPerNodes(innerNode, :); % Trovo gli indici dei triangoli che hanno innerNode come Vertice
outerVertexes = indSMatrix(outerNode,:); % Trovo che vertice è (primo secondo terzo) outerNode per ogni triangolo di 
                                                % outerTriangles
innerVertexes = indSMatrix(innerNode,:); % Stessa cosa per il nodo interno


for i = 1 : maxNumTriangles % Ciclo sui triangoli esterni
    
    currentOuterTriangle = outerTriangles(i); % Triangolo esterno corrente

    if (currentOuterTriangle == 0) % se è zero vuol dire che li abbiamo già fatti tutti
        break;
    end
    
     % Voglio individuare che vertice è outerNode rispetto a
        % currentOuterTriangle e che vertice è innerNode rispetto a
        % CurrentInnerTriangle in modo da estrarre  le "vele" giuste per il
        % calcolo del nucleo. La "vela"giusta sul triangolo corrente è
        % quella che fa 1 sul nodo corrente, quindi rispetto al triangolo,
        % quella che fa 1 sul vertice corrente

    currentOuterVertex = outerVertexes(i); % il punto geometrico a cui ci riferiamo è sempre lo stesso, ma cambia in quanto
    % vertice dei triangoli che toccano tale punto

    for j = 1 : maxNumTriangles %Ciclo sui triangoli interni
        currentInnerTriangle = innerTriangles(j); % Triangolo interno corrente

        if (currentInnerTriangle == 0) % Se è zero vuol dire che li abbiamo già fatti tutti
        break;
        end
        
        

        currentInnerVertex = innerVertexes(j);
        %per mappare i nodi di Gauss-Hammer sui triangoli mi serviranno i
        %vertici dei due triangoli per cui li estraggo qui in due matrici
        %3x3 

       


        currentOuterVerts = zeros(3,3); %Nella prima riga metto le 3 coordinate del primo vertice e così via

        currentOuterVerts(1,:) = domainMesh.coordinates(domainMesh.triangles(currentOuterTriangle,1),:);
        currentOuterVerts(2,:) = domainMesh.coordinates(domainMesh.triangles(currentOuterTriangle,2),:);
        currentOuterVerts(3,:) = domainMesh.coordinates(domainMesh.triangles(currentOuterTriangle,3),:);

        currentInnerVerts = zeros(3,3); % Initialize matrix for inner triangle vertices

        currentInnerVerts(1,:) = domainMesh.coordinates(domainMesh.triangles(currentInnerTriangle,1),:);
        currentInnerVerts(2,:) = domainMesh.coordinates(domainMesh.triangles(currentInnerTriangle,2),:);
        currentInnerVerts(3,:) = domainMesh.coordinates(domainMesh.triangles(currentInnerTriangle,3),:);

        % Per il nucleo mi servono anche le normali dei triangoli esterno
        % ed intrno correnti

        
        currentOuterNormal = domainMesh.normal(currentOuterTriangle,:); % Normale del triangolo esterno

        
        currentInnerNormal = domainMesh.normal(currentInnerTriangle,:); % Normale del triangolo interno

        
        % Ora estraggo i valori necessari per calcolare la funzione di base corrente ed il vettore
        % associato al suo rotore
        
        currentOuterTriangleMatrix = constData{currentOuterTriangle}.matCoeff;
        currentInnerTriangleMatrix = constData{currentInnerTriangle}.matCoeff;
        
        currentOuterTriangleCoeffs = constData{currentOuterTriangle}.vetCoeff; %Il pezzo della "vela" che dipende solo dal triangolo e non dal punto o dal vertice
        currentInnerTriangleCoeffs = constData{currentInnerTriangle}.vetCoeff; %Il pezzo della "vela" che dipende solo dal triangolo (interno)

        currentOuterVVector = cross(currentOuterNormal,currentOuterTriangleMatrix(currentOuterVertex,:)); % è \bm{V_{\tilde{\alpha}}^{\tilde{s}}}
        currentInnerVVector = cross(currentInnerNormal, currentInnerTriangleMatrix(currentInnerVertex,:)); % è \bm{V_{\alpha}^{s}}
        

        if (currentOuterTriangle == currentInnerTriangle)
            singularSubBlock = computeSingBlockW(pbParam, currentOuterVertex, currentInnerVertex, currentOuterNormal, constData{currentOuterTriangle}, timeInstant, quadData.methodSpecs, delta, epsilon);
            OutputMatrix = OutputMatrix + singularSubBlock; 
            continue; % Salta tutto il resto del ciclo e passa al prossimo 'j'
        end



        % Ora dovrei avere tutto e posso cominciare la quadratura facendo
        % il doppio ciclo sui nodi di gauss (su GPU il ciclo esterno sarà gestito dai Thread nel blocco corrente, quello interno invece sarà sempre un for)

        for innerGaussHammerIndex = 1 : length(quadData.INTn)
            %devo calcolare il peso e le coordinate del dodo corrente
            currentInnerGaussWeight = quadData.INTw(mod(innerGaussHammerIndex,3)+1)*domainMesh.area(currentInnerTriangle); % Peso del nodo corrente sul triangolo corrente
            %INtw e INTn li tirerò fuori da quad data
            temporaryStandardNode = quadData.INTn(innerGaussHammerIndex,:); % Nodo su triangolo di riferimento

            currentInnerGaussNode = temporaryStandardNode*currentInnerVerts; % Mappo il nodo sul triangolo interno corrente

            currentVXi = VFunction(currentInnerTriangleCoeffs, currentInnerTriangleMatrix, currentInnerVertex, currentInnerGaussNode); % Funzione che creerò per valutare V_{\alpha}^{s}(\xi)

            for outerGaussHammerIndex = 1 : length(quadData.EXTn)
                currentOuterGaussWeight = quadData.EXTw*domainMesh.area(currentOuterTriangle); % Peso del noto esterno

                temporaryStandardNode = quadData.EXTn(outerGaussHammerIndex,:); %Nodo sul triangolo di riferimento (esterno)

                currentOuterGaussNode = temporaryStandardNode*currentOuterVerts; % Mappo il nodo sul triangolo Esterno Corrente

                currentVTildeX = VTildeFunction(currentOuterTriangleCoeffs, currentOuterTriangleMatrix, currentOuterVertex, currentOuterGaussNode); % Funzione che valuta V_{\tilde{\alpha}}^{\tilde{s}}(x)
                    

                current_rVector = currentOuterGaussNode-currentInnerGaussNode; %vettore differenza r = x-\xi

                current_rNorm = norm(current_rVector,2); % norma di r
                
                
                % Ciclo sugli istanti temporali

                for k = 1 : 4
                    currentTime = (timeInstant+timeCoeffs(k))*deltaT;
                    if (currentTime<=0)
                        continue;
                    end
                    currentCoeff = kernelCoeffs(k);
                    kernel = kernelCalc(currentTime, currentVTildeX,currentVXi,currentInnerVVector,currentOuterVVector,...
                        current_rVector, current_rNorm, currentInnerNormal, currentOuterNormal, mu, lambda, rho, cS, cP, delta, epsilon);

                    OutputMatrix = OutputMatrix + currentInnerGaussWeight*currentOuterGaussWeight*(currentCoeff/(deltaT^2))*kernel;
                        
                end

            end
        end

    end

end

end



function VXi = VFunction(triangleCoeffs, triangleMatrix, vertexInd, xi)
%% Questa funzione serve per calcolare V_{\alpha}^{s}(\xi)
eVec = [0,0,0];
eVec(vertexInd) = 1;
VXi = dot((triangleCoeffs + triangleMatrix*xi'), eVec); %traspongo xi perché viene passato come riga
end

function VTildeX = VTildeFunction(triangleCoeffs, triangleMatrix, vertexInd, x)
%% Questa funzione serve per clcolare V_{\tilde{\alpha}}^{\tilde{s}}(x)
eVec = [0,0,0];
eVec(vertexInd) = 1;
VTildeX = dot((triangleCoeffs +triangleMatrix*x'), eVec); %traspongo x perché viene passato come riga
end

function kernelTot = kernelCalc(t, VTildeX, VXi, VAlphaS, VTildeAlphaS, rVector, r, n, v, mu, lambda, rho, cS, cP, delta, epsilon)
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

kernelTDelta1 = calcKernelTDelta1(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta);
kernelTDelta2 = calcKernelTDelta2(i, k, t, VTildeX, VXi, r, n, v, mu, lambda, rho, cS, cP, delta);
kernelTH = calcKernelTH(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta);

kernelT = kernelTDelta1 + kernelTDelta2 + kernelTH;
end


function kernelTDelta1 = calcKernelTDelta1(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta)
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
H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IDeltaS2 = H_S/(cS*cS*cS*cS); 
IDeltaP2 = H_P/(cP*cP*cP*cP);

Nv = dot(n,v);

constTerm = -(lambda*mu + 2*mu*mu)/(4*pi*rho*(lambda+mu));

kernelTDelta2 = VTildeX*VXi*constTerm*((lambda*n(k)*v(i) + mu*n(i)*v(k) + mu*delta(i,k)*(Nv))/(r))*(IDeltaS2-IDeltaP2);
end

function kernelTH = calcKernelTH(i, k, t, VTildeX, VXi, rVector, r, n, v, mu, lambda, cS, cP, delta)

H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IHS = H_S*(t/2);
IHP = H_P*(t/2);

Rn = dot(rVector,n);
Rv = dot(rVector,v);
Nv = dot(n,v);

constTerm = mu/(2*pi*(lambda+mu));

kernelTH = VTildeX*VXi*constTerm*(3*lambda*((rVector(k)*v(i)*Rn + rVector(i)*n(k)*Rv)/(r*r*r*r*r)) - 2*lambda*((n(k)*v(i))/(r*r*r))...
               +3*mu*((rVector(k)*n(i)*Rv + rVector(i)*v(k)*Rn + delta(i,k)*Rv*Rn + rVector(i)*rVector(k)*Nv)/(2*r*r*r*r*r)) - mu*((n(i)*v(k) + delta(i,k)*Nv)/(r*r*r)))*(IHS-IHP);

end


function kernelR = kernelRCalc(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon)

kernelRDelta = calcKernelRDelta(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon);
kernelRH = calcKernelRH(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon);

kernelR = kernelRH + kernelRDelta;

end

function kernelRH = calcKernelRH(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon)
kernelRH = 0;
H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IHRGS = H_S*(((t-(r/cS))^4)/(24) - r*((t-(r/cS))^3)/(6*cS));
IHRGP = H_P*(((t-(r/cP))^4)/(24) - r*((t-(r/cP))^3)/(6*cP));

constTerm = -((mu*mu)/(4*pi*rho*(lambda+mu)));

for m = 1 : 3
    for q = 1 : 3
        for b = 1 : 3
            for g = 1 : 3
                for e = 1 :3
                    for s = 1 : 3

                        kernelRH = kernelRH + VTildeAlphaS(m)*VAlphaS(q)*(epsilon(i,b,e)*epsilon(k,g,s)...
                            *(2*lambda*delta(e,m)*delta(s,q) + (lambda+2*mu)*(delta(e,q)*delta(m,s) + delta(e,s)*delta(m,q)))...
                            *(3*((rVector(b)*rVector(g))/(r*r*r*r*r)) - delta(b,g)/(r*r*r)));

                    end
                end
            end
        end
    end
end
kernelRH = constTerm*kernelRH*(IHRGS-IHRGP);
end

function kernelRDelta = calcKernelRDelta(i, k, t, VAlphaS, VTildeAlphaS, rVector, r, mu, lambda, rho, cS, cP, delta, epsilon)
kernelRDelta = 0;
H_P = (t - (r/cP) > 0); % queste sono le heavyside
H_S = (t - (r/cS) > 0);

IDeltaRGS = H_S*(((t-(r/cS))^2)/(2*cS*cS));
IDeltaRGP = H_P*(((t-(r/cP))^2)/(2*cP*cP));

constTerm = -((mu*mu)/(4*pi*rho*(lambda+mu)));

for m = 1 : 3
    for q = 1 : 3
        for b = 1 : 3
            for g = 1 : 3
                for e = 1 :3
                    for s = 1 : 3

                        kernelRDelta = kernelRDelta + VTildeAlphaS(m)*VAlphaS(q)*(epsilon(i,b,e)*epsilon(k,g,s)...
                            *(2*lambda*delta(e,m)*delta(s,q) + (lambda+2*mu)*(delta(e,q)*delta(m,s) + delta(e,s)*delta(m,q)))...
                            *((rVector(b)*rVector(g))/(r*r*r)));

                    end
                end
            end
        end
    end
end
kernelRDelta = constTerm*kernelRDelta*(IDeltaRGS-IDeltaRGP);

end

function singularSubBlock = computeSingBlockW(pbParam, outerVertex, innerVertex, triangleNormal, currentConstData, timeInstant, methodSpecs, delta, epsilon)

singularSubBlock = zeros(3,3);

mu = pbParam.mu;
lambda = pbParam.lambda;
rho = pbParam.rho;
deltaT = pbParam.Tfin/pbParam.nT;
cS = pbParam.velS;
cP = pbParam.velP;

kernelCoeffs = [1,-3,3,1]; % I coefficienti C_\eta
timeCoeffs = [-2,-1,0,1]; % Gli \eta che individuano i vari istanti temporali
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

            innerIntegral = zeros(3,3);

            for indChild = 1 : 3 %ciclo sui triangoli figli interni
                rMin = max(cS * minDiffTemp, 0);
                rInt = cS * currentTime;
                rExt = cP * currentTime;

                [G2DCnodes, G2DCweights] = eebem.utility.quadratureRules.generateFinalG2Dnodes(currentConstData.childVerts{indEXTn, indChild}, rMin, rInt, rExt, methodSpecs.numSNGLR);
                
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




