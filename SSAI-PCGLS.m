%% load matrix
clear all;
% well-conditioned
%load('image_interp.mat');
%load('Hardesty2.mat');
%load('Hardesty3.mat');
%load('ch8-8-b2.mat');
%load('ch8-8-b3.mat');
%load('ch8-8-b4.mat');
%load('ch8-8-b5.mat');
%load('LargeRegFile.mat');
%load('Rucci1.mat');
load('sls.mat');
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
x_0=spalloc(nmax,1,nmax);
%% normalize A
C=diag(sparse(1./sqrt(diag(A2))));
A2=C*tril(A2,-1)*C;
A2=A2+A2'+speye(nmax);
tol=1e-8*norm(b);
if (isnan(tol))
    exit;
end
%% sparse inverse
lfil=ceil(nnz(A2)/nmax);
tic
M=entire_r_sparse_inverse(A2,nmax,lfil);
M=(M+M')/2;
toc;
%% PCGLS
tic;
[x,iter,flag]=spcgls(A,b,tol,nmax,M);
toc;
disp(iter)