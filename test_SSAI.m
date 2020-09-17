%this file just shows how to call the function
clear all;
load('cvxbqp1.mat');
A=Problem.A;
n=length(A);
b=sparse(A*(1:n)'/n);
x_0=spalloc(n,1,n);
tol=1e-8;
[x,err,res,iter]=SSAI_PCG(A,b,x_0,tol);