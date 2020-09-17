%% load matrix
clear all;
load('ch8-8-b5.mat');
%% problem setup
A=Problem.A;
onnz=nnz(A);
row=length(A);
%% solution setup
b=ones(row,1);
%%
A2=A'*A;
nnnz=nnz(A);
tol2=1e-2;
nmax=length(A2);
%% normalize A
C=diag(sparse(1./sqrt(diag(A2))));
A2=C*tril(A2,-1)*C;
A2=A2+A2'+speye(nmax);
tol=1e-8;
%% sparse inverse
lfil=ceil(nnz(A2)/nmax);
tic
M=entire_r_sparse_inverse(A2,nmax,lfil);
M=(M+M')/2;
toc;
%% SPCGLS
tic;
[x,iter,flag]=spcgls(A*C,b,tol*norm(b),nmax,M);
toc;
disp(iter)

