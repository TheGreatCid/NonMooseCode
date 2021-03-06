%Kl Expansion
%David Torres
clc;clear;close all;
%% Common Values
correlation = .020; %mm
standdev = .05; % 5%
calc = 1;
%% Importing Mesh

fid = fopen('thickwall.msh');

hdrRows = 11;
hdrData = textscan(fid,'%s', hdrRows,'Delimiter','\n');
matData = textscan(fid,'%f%f%f%f','Delimiter',' ','CollectOutput',true);
hdrData2 = textscan(fid,'%s', 121, 'Delimiter','\n');
matData2 = textscan(fid,'%f%f%f%f%f%f%f%f','Delimiter',' ','CollectOutput',true);

fclose(fid);
nodes = matData{1}(:,2:3);
el = matData2{1}(:,6:8)';



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
xi1 = [1/6; 2/3; 1/6]; %Integration points
xi2 = [1/6; 1/6; 2/3];
w = [1/6; 1/6; 1/6]; %Weighting
lam = 1-xi1-xi2; %phi1 prewritten to save space
C = zeros(cN,cN);

Nij = zeros(cN,cN);

Q = 3;
if calc == 1
for m = 1:rE %Outer element loop
    CeOuter = 0; %Reset variable - Outer set of summations\
    Jm = [nodes(el(1,m),1)-nodes(el(3,m),1) nodes(el(1,m),2)-nodes(el(3,m),2);...
        nodes(el(2,m),1)-nodes(el(3,m),1) nodes(el(2,m),2)-nodes(el(3,m),2) ]; %jacobian
    idx2 = [el(1,m), el(2,m), el(3,m)];
    %Ne = 0;

    for p = 1:Q %Outer qp loop        
        
        Np = [lam(p) xi1(p) xi2(p)];... %TRI3 shape functions at current integration points
            
        x2 = [ Np(1)*nodes(el(1,m),1) + Np(2)*nodes(el(2,m),1) + Np(3)*nodes(el(3,m),1);... %2nd coordinate
            Np(1)*nodes(el(1,m),2) + Np(2)*nodes(el(2,m),2) + Np(3)*nodes(el(3,m),2)];
        
        
        for n = 1:rE %Inner Element loop
            
            CeInn = 0; %Reset variable - Inner set of summations
            Jn = [nodes(el(1,n),1)-nodes(el(3,n),1) nodes(el(1,n),2)-nodes(el(3,n),2);nodes(el(2,n),1)-nodes(el(3,n),1) nodes(el(2,n),2)-nodes(el(3,n),2) ]; %jacobian
            
            for q = 1:Q %Inner qp loop
                
                Nq = [lam(q) xi1(q)  xi2(q)]; %TRI3 shape functions at current integration points
                
                x1 =[ Nq(1)*nodes(el(1,n),1) + Nq(2)*nodes(el(2,n),1) + Nq(3)*nodes(el(3,n),1);... %1st coordinate
                    Nq(1)*nodes(el(1,n),2) + Nq(2)*nodes(el(2,n),2) + Nq(3)*nodes(el(3,n),2)];
                
                
                Cs = standdev^2*exp(-(sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2))/correlation);
                
                
                CeInn = CeInn + Cs*Nq'*abs(det(Jn))*w(q); %Calculate inner summation 3x3
                
                
            end
            
            %Assemble CeInn Value into global matrix
            idx = [el(1,n), el(2,n), el(3,n)]; %Assemble over element n
            C(idx,idx2) = C(idx,idx2) + CeInn*Np*abs(det(Jm))*w(p); %Put values where nodal DoFs are located
            
        end        
        
        
       % Ne = Ne + Np*Np'*det(Jm)*w(p); %Calc N matrix
        
        
    end
    
   % Nij(idx2,idx2) = Nij(idx2,idx2) + Ne; %Put values where nodal DoFs are located
    %Integrating about single variable but still to loop about both sets of
    %elements 
    if rem(m,100)==0
        clc;
        fprintf('%1.1f%%\n',m/rE*100)
    end
    
end

%% Calc N

for m = 1:rE
    Jm = [nodes(el(1,m),1)-nodes(el(3,m),1) nodes(el(1,m),2)-nodes(el(3,m),2);nodes(el(2,m),1)-nodes(el(3,m),1) nodes(el(2,m),2)-nodes(el(3,m),2) ]; %jacobian

   % for n = 1:rE
       % Ne = 0;
        %Jn = [nodes(el(1,n),2)-nodes(el(3,n),2) nodes(el(1,n),2)-nodes(el(3,n),2);nodes(el(2,n),1)-nodes(el(3,n),1) nodes(el(2,n),1)-nodes(el(3,m),1) ]; %jacobian
        Ne = 0;
        for p=1:Q
            
            Np = [lam(p) xi1(p) xi2(p)];
            
            Ne = Ne + Np'*Np*abs(det(Jm))*w(p);
            
        end
        idx = [el(1,m), el(2,m), el(3,m)]; %Assemble over element n
        idx2 = [el(1,m), el(2,m), el(3,m)];

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
save('EigVal.mat','D');
save('EigVec.mat','V');
else
    
    D = load('EigVal.mat');
    V = load('EigVec.mat');
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
trisurf(tri,nodes(:,1),nodes(:,2),f)
colorbar








