%Kl Expansion
%David Torres
clc;clear;close all;
%% Common Values
correlation = .020; %m
standdev = 1; % 5%
calc = 0;
%% Importing Mesh

fid = fopen('KLquarter.msh');

hdrRows = 13;
hdrData = textscan(fid,'%s', hdrRows,'Delimiter','\n');
matData = textscan(fid,'%f%f%f%f','Delimiter',' ','CollectOutput',true);
hdrData2 = textscan(fid,'%s', 123, 'Delimiter','\n');
matData2 = textscan(fid,'%f%f%f%f%f%f%f%f%f','Delimiter',' ','CollectOutput',true);

fclose(fid);
nodes = matData{1}(:,2:3);
el = matData2{1}(:,6:end)';



%% Create Shape Functions
%local shape function -
%2-D linear triangles

%x = phi1(xi,nu)x1 + phi2(xi,nu)x2 + phi3(xi,nu)x3
%y = phi1(xi,nu)y2 + phi2(xi,nu)y2 + phi3(xi,nu)y3
%phi1(x,y) = a1 + b1*xi + c1*nu

%phii = ai+bix+ciy
[NumNodes, ~] = size(nodes);
[~,NumEl] = size(el);
%N = zeros(cN);

%fprintf('\n Done! \n')

%% Calc Values

%xi1 = [1/sqrt(3); -1/sqrt(3)]; %Integration points
%xi2 = [1/sqrt(3); -1/sqrt(3)];

xi = [1/sqrt(3), -1/sqrt(3), -1/sqrt(3), 1/sqrt(3);
      1/sqrt(3), 1/sqrt(3), -1/sqrt(3), -1/sqrt(3)];

w = [1 1 1 1];
   %Weighting

C = zeros(NumNodes,NumNodes);
Nij = zeros(NumNodes,NumNodes);
Q = 4;
if calc == 1

  
 %Parametric shape functions
for m = 1:NumEl %Outer element loop        
        idxm = el(:,m);
        xm = nodes(idxm,:); %coordinatates for mth element
             
        %Ne = 0;
        
        
    

         Ne = 0;
        for p = 1:Q %Outer qp loop 1
           
            Np = shape(xi(:,p));
            coordp = (Np*xm)';
            
            Gp = gradShape(xi(:,p)); %gradient of shape functions
            Jp = Gp*xm; %Jacobian
           
            
          
            
            for n = 1:NumEl %Inner Element loop
                
                CeInn = 0; %Reset variable - Inner set of summations
                idxn = el(:,n);
                xn = nodes(idxn,:);
                
                for q = 1:Q %Inner qp loop
                    
                    
                    Nq = shape(xi(:,q));
                    coordq = (Nq*xn)';
                    
                    Gq = gradShape(xi(:,q)); %gradient of shape functions
                    Jq = Gq*xn; %Jacobian
                    
                    
                    
                    
                    tau = norm(coordp-coordq);
                    Cs = standdev^2*exp(-tau^2/correlation^2);
                    
                    CeInn = CeInn + Cs*Nq*det(Jq)*w(q); %Calculate inner summation                     
                end

                %Assemble CeInn Value into global matrix
                C(idxm,idxn) = C(idxm,idxn) + Np'*CeInn*det(Jp)*w(p); %Put values where nodal DoFs are located
           
            end
            
             Ne = Ne + Np'*Np*det(Jp)*w(p);

        end
            Nij(idxm,idxm) = Nij(idxm,idxm) + Ne;
        
        % Ne = Ne + Np*Np'*det(Jm)*w(p); %Calc N matrix
        
        % Nij(idx2,idx2) = Nij(idx2,idx2) + Ne; %Put values where nodal DoFs are located
        %Integrating about single variable but still to loop about both sets of
        %elements
        if rem(m,10)==0
            clc;
            fprintf('%1.1f%%\n',m/NumEl*100)
        end
        
end
        
        

        
 
    
    
    
    
    spy(Nij)
    figure
    spy(C)
    
    
    %% Create EigenValues
    
    [V,D] = eig(C,Nij); %Takes a long time, so only run when I know above
    %calcs are right
    %D consists of generalized eigenvalues
    %V are the corresponding eigenvectors
    s = sum(D);
    s50 = s(1:50);
    figure
    plot(1:50,s50,'.')
    xlim([1,50])
    %%
    save('EigValQuart.mat','D');
    save('EigVecQuart.mat','V');
else
    
    D = load('EigValQuart.mat');
    V = load('EigVecQuart.mat');
    D = D.D;
    V = V.V;
end

%%
[D,idx] = sort(diag(D),'descend');

D2 = sqrt(D(1:50)); %squrt of first 50 eigen values

V2 = V(:,idx); %First 50 eigenvectors
V2 = V2(:,1:50);

eta = randn(50,1);
Xi = V2*(D2.*eta); %First 50 eigenvectors * First 50 eigenvectors * Normalized distibution
%normcdf
%conver to gamme dist, gampdf
E = 210e9; %Mean young's modulus

A = (E^2+standdev*2)/standdev;
B = E*(A-1);

Phi = normcdf(Xi); %Calc phi
P = gaminv(Phi,A,B); %Calc p
%% Contour;

dt = delaunayTriangulation(nodes(:,1),nodes(:,2));
tri = dt.ConnectivityList;
figure
%trisurf(tri,nodes(:,1),nodes(:,2),f)
%trisurf(tri,nodes(:,1),nodes(:,2),real(f))
colorbar
for m = 1:NumEl
    idxm = el(:,m);
    xm = nodes(idxm,:);
    patch(xm(:,1)',xm(:,2)',Xi(idxm));
    hold on
end
title('Patch')
figure
for m = 1:NumEl
    idxm = el(:,m);
    xm = nodes(idxm,:);
    patch(xm(:,1)',xm(:,2)',P(idxm));
    hold on
end
colorbar
title('Patch')
%% Get values from FE mesh
fid = fopen('quarter.msh');

hdrRows = 13;
FEhdrData = textscan(fid,'%s', hdrRows,'Delimiter','\n');
FEmatData = textscan(fid,'%f%f%f%f','Delimiter',' ','CollectOutput',true);
FEhdrData2 = textscan(fid,'%s', 1198, 'Delimiter','\n');
FEmatData2 = textscan(fid,'%f%f%f%f%f%f%f%f%f','Delimiter',' ','CollectOutput',true);

fclose(fid);
FEnodes = matData{1}(:,2:3);
FEel = matData2{1}(:,6:end)';

[FEcN, FErN] = size(FEnodes);
[FEcE,FErE] = size(FEel);

% %% Transfer values to FE mesh
% 
% C_FE = zeros(FEcN);
% %xis(1,:) = [0,0,0,0];
% %xis(2,:) = [0,0,0,0];
% m = 1;
% fun = @(x)sums(x(1,:),x(2,:),m,FEel,el,nodes,FEnodes);
% for j = 1:FErE
%    
%     
%     idx = [el(1,j), el(2,j), el(3,j), el(4,j)];
%     %F(:,:) = sums(xis,m,FEel,el,nodes,FEnodes);
%     xis = fsolve(fun,[0,0,0,0;0,0,0,0]);
%     nodepara(m,:) = [xis(1,:);xis(2,:)];
%     
%     
% end










