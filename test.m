%this file just shows how to call the function
A=sparse([2,1,1;1,3,1;1,1,5]);
b=sparse(3,1,1);
n=length(A);
x_0=spalloc(n,1,n);
tol=1e-8;
[x,err,res,iter]=SSAI(A,b,x_0,tol);