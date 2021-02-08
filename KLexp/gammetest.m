clc;clear;

%mean 15
%variance .1
A = 1/.1;
B = 15/A;
phi = randn(2000,1);

Phi = normcdf(phi); %Calc phi
P = gaminv(Phi,A,B); %Calc p

