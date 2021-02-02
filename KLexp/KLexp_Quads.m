%Kl Expansion
%David Torres
clc;clear;close all;
%% Common Values
correlation = .020; %m
standdev = .05; % 5%
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
[cN, rN] = size(nodes);
[cE,rE] = size(el);
%N = zeros(cN);

%fprintf('\n Done! \n')

%% Calc Values
n = 3; %nodes per element
xi1 = [1/sqrt(3); -1/sqrt(3)]; %Integration points
xi2 = [1/sqrt(3); -1/sqrt(3)];
w = [1; 1]; %Weighting
lam = 1-xi1-xi2; %phi1 prewritten to save space
%save C.mat
%D = zeros(cN,cN);
C = zeros(cN,cN);
Nij = zeros(cN,cN);
Q = 2;
if calc == 1

  
 %Parametric shape functions
    for m = 1:rE %Outer element loop
        CeOuter = 0; %Reset variable - Outer set of summations\
        
        xem = nodes(el(:,m),1); %xe corrdinatates for mth element
        yem = nodes(el(:,m),2); %ye coordinates for mth element
        idx2 = [el(1,m), el(2,m), el(3,m) el(4,m)];
        %Ne = 0;
        
        
        
        
        for p = 1:Q %Outer qp loop 1
           
               
        
        
        
        for k = 1:Q %Outer qp loop 2
            x1 = xem(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + xem(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + xem(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + xem(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
            y1 = yem(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + yem(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + yem(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + yem(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
            G = [xi2(k)-1 1-xi2(k) xi2(k)+1 -xi1(k)-1;...
                xi1(p)-1 -xi1(p)-1 xi1(p)+1 1-xi1(p)]; %gradient of shape functions
            Jm = G*[xem,yem]; %Jacobian
            Npk = [(1/4)*(1-xi1(p))*(1-xi2(k)),(1/4)*(1+xi1(p))*(1-xi2(k)),(1/4)*(1+xi1(p))*(1+xi2(k)),(1/4)*(1-xi1(p))*(1+xi2(k))];
            
          
            
            for n = 1:rE %Inner Element loop
        
                CeInn = 0; %Reset variable - Inner set of summations
                
                for q = 1:Q %Inner qp loop
                    for l = 1:Q %Inner qp loop 
                       Nql  = [(1/4)*(1-xi1(q))*(1-xi2(l)),(1/4)*(1+xi1(q))*(1-xi2(l)),(1/4)*(1+xi1(q))*(1+xi2(l)),(1/4)*(1-xi1(q))*(1+xi2(l))];
                         xen = nodes(el(:,n),1); %xe corrdinatates for nth element
                         yen = nodes(el(:,n),2); %ye coordinates for nth element
                         
                        G = [xi2(l)-1 1-xi2(l) xi2(l)+1 -xi1(l)-1;...
                                 xi1(q)-1 -xi1(q)-1 xi1(q)+1 1-xi1(q)];
                        Jn = G*[xen,yen];
                                     
                        x2 = xen(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + xen(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + xen(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + xen(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
                        y2 = yen(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + yen(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + yen(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + yen(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
                        Cs = standdev^2*exp(-(sqrt((x1-x2)^2+(y1-y2)^2))/correlation);
                        
                        CeInn = CeInn + Cs*Nql*abs(det(Jn))*w(q)*w(l); %Calculate inner summation 3x3
                        
                    end
                end
                
                %Assemble CeInn Value into global matrix
                idx = [el(1,n), el(2,n), el(3,n) el(3,n)]; %Assemble over element n
                C(idx,idx2) = C(idx,idx2) + CeInn'*Npk*abs(det(Jm))*w(p)*w(k)'; %Put values where nodal DoFs are located
           
            end
        end
        
        % Ne = Ne + Np*Np'*det(Jm)*w(p); %Calc N matrix
        
        
        end
        
        % Nij(idx2,idx2) = Nij(idx2,idx2) + Ne; %Put values where nodal DoFs are located
        %Integrating about single variable but still to loop about both sets of
        %elements
        if rem(m,10)==0
            clc;
            fprintf('%1.1f%%\n',m/rE*100)
        end
        
    end
    
    %% Calc N
    
    for m = 1:rE
        x1 = xem(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + xem(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + xem(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + xem(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
        y1 = yem(1)*(1/4)*(1-xi1(p))*(1-xi2(k)) + yem(2)*(1/4)*(1+xi1(p))*(1-xi2(k)) + yem(3)*(1/4)*(1+xi1(p))*(1+xi2(k)) + yem(4)*(1/4)*(1-xi1(p))*(1+xi2(k));
        G = [xi2(k)-1 1-xi2(k) xi2(k)+1 -xi1(k)-1;...
             xi1(p)-1 -xi1(p)-1 xi1(p)+1 1-xi1(p)]; %gradient of shape functions
        Jm = G*[xem,yem]; %Jacobian
        
        % for n = 1:rE
        % Ne = 0;
        %Jn = [nodes(el(1,n),2)-nodes(el(3,n),2) nodes(el(1,n),2)-nodes(el(3,n),2);nodes(el(2,n),1)-nodes(el(3,n),1) nodes(el(2,n),1)-nodes(el(3,m),1) ]; %jacobian
        Ne = 0;
        for p=1:Q
            for k = 1:k
            
               Npk = [(1/4)*(1-xi1(p))*(1-xi2(k)),(1/4)*(1+xi1(p))*(1-xi2(k)),(1/4)*(1+xi1(p))*(1+xi2(k)),(1/4)*(1-xi1(p))*(1+xi2(k))];
            
             Ne = Ne + Npk'*Npk*abs(det(Jm))*w(p)*w(k);
            end
        end
        idx = [el(1,m), el(2,m), el(3,m), el(4,m)]; %Assemble over element n
        idx2 = [el(1,m), el(2,m), el(3,m), el(4,m)];
        
        Nij(idx,idx2) = Nij(idx,idx2) + Ne; %Put values where nodal DoFs are located
        %end
        if rem(m,100)==0
            clc;
            fprintf('%1.1f%%\n - N',m/rE*100)
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

%% Use eigenvalues
% syms x y
% Create spacial shape functions
%
% for i = 1:rE
%     x1 = nodes(el(1,i),1);
%     x2 = nodes(el(2,i),1);
%     x3 = nodes(el(3,i),1);
%     y1 = nodes(el(1,i),2);
%     y2 = nodes(el(2,i),2);
%     y3 = nodes(el(3,i),2);
%    count = count+1;
%
% %     N(1,i) = (1/(2*A(i)))*(x2*y3-x3*y2+(y2-y3)*x + (x3-x2)*y);
% %         N(2,i) = (1/(2*A(i)))*(x3*y1-x1*y3+(y3-y1)*x + (x1-x3)*y);
% %         N(3,i) = (1/(2*A(i)))*(x1*y2-x2*y1+(y1-y2)*x + (x2-x1)*y);
%      ind1 = el(1,i);
%      ind2 = el(2,i);
%      ind3 = el(3,i);
%
%
%     A(i) = abs(det([x2-x1 y2-y1;x3-x1 y3-y1])*1/2); % Calculating each element area for later use
%     for j = 1:3 %Create mult 1/(2Ae) as 1st index in 3D matrix
%        N(ind1,i,1) = 1/(2*A(i));
%        N(ind2,i,1) = 1/(2*A(i));
%        N(ind3,i,1) = 1/(2*A(i));
%     end
%
%     for j = 1:3 %Create 1st constant value as 2nd value in 3D matrix
%         N(ind1,i,2) = x2*y3;
%         N(ind2,i,2) = x3*y1;
%         N(ind3,i,2) = x1*y2;
%     end
%
%     for j = 1:3 %Create const mult by x as 3rd value in 3D matrix
%        N(ind1,i,3) = y2-y3;
%        N(ind2,i,3) = y3-y1;
%        N(ind3,i,3) = y1-y2;
%     end
%
%     for j = 1:3 %Create const mult by y as 4th value in 3D matrix
%        N(1,i,4) = x3-x2;
%        N(2,i,4) = x1-x3;
%        N(3,i,4) = x2-x1;
%     end
%
%
%
%     if rem(i,100) == 0
%         clc;
%         fprintf('%1.1f%% - Making Spatial Shape Functions',i/rE*100)
%     end
% end


%Calculate eigen functions


%don't spatial shape function equal 1 at nodes?

% f = zeros(1,cN);
% for k = 1:cN
%     for j = 1:cN
%         for i = 1:3
%             %ind = find(el==j)
%             f(k) = f(k) + V(j,k);
%         end
%     end
% end
% %S(x,theta)
% %normrnd(mu,sigma)
% %use M=50
% M = 50;
% Ms = 210;
%
% eig = s50;
%
% for i = 1:cN
%        theta = normrnd(100,50);
%     for k = 1:50
%
%
%         s(i) = s(i) + real(sqrt(D(k))*f(k)*theta);
%     end
%
% end
% s = real(s);

D2 = sqrt(D(1:50, 1:50));
V2 = V(:,1:50);

xi = normrnd(zeros(1,50),standdev*100)*210e9;
XI = diag(xi);
DV = V2*D2*XI;
f = sum(DV,2)+210e9*ones(cN,1);
figure
%contour(nodes(:,1),nodes(:,2),real(s))
%
% xv = min(nodes(:,1)):max(nodes(:,1));
% yv = min(nodes(:,2)):max(nodes(:,2));
% %[xv,yv] = meshgrid(nodes(:,1),nodes(:,2));
% [X,Y] = ndgrid(xv(1,:),yv(1,:));
% Z = griddata(nodes(:,1),nodes(:,2),real(s),X,Y);
% figure
% contour(X,Y,Z,100)
%% Contour

dt = delaunayTriangulation(nodes(:,1),nodes(:,2));
tri = dt.ConnectivityList;
figure
%trisurf(tri,nodes(:,1),nodes(:,2),f)
trisurf(tri,nodes(:,1),nodes(:,2),real(f))
colorbar


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

%% Transfer values to FE mesh

C_FE = zeros(FEcN);
%xis(1,:) = [0,0,0,0];
%xis(2,:) = [0,0,0,0];
m = 1;
fun = @(x)sums(x(1,:),x(2,:),m,FEel,el,nodes,FEnodes);
for m = 1:FErE
   
    
    idx = [el(1,m), el(2,m), el(3,m), el(4,m)];
    %F(:,:) = sums(xis,m,FEel,el,nodes,FEnodes);
    xis = fsolve(fun,[0,0,0,0;0,0,0,0]);
    C_FE(idx,idx) = C_FE(idx,idx) + [xis(1,:),xis(2,:)];
    
    
end










