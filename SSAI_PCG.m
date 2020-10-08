function [x_0,err,res,iter]=SSAI_PCG(A,b,x_0,tol)
% Solves a linear HPD system Ax=b, with initial guess x_0
% b and x_0 are column vectors
% tol is the relative tolerance to inital guess to solve for
%x_0 is overwritten with the solution
% err is the residual norm
% res is the number of restarts, iter is the number of iterations
tol2=1e-2;
n=length(A);
%% normalize A
C=diag(sparse(1./sqrt(diag(A))));
A=C*tril(A,-1)*C;
A=A+A'+speye(n);
%% solution setup
tol=tol*norm(b);
%% sparse inverse
lfil=ceil(nnz(A)/n);
tic
M=entire_r_sparse_inverse(A,n,lfil);
M=(M+M')/2;
toc;
%% PCG
delta=10;
waste_iter=0;
tic;
flag=1;
count=0;
while(flag==1)
    dx=spalloc(n,1,n);
    r=b-A*x_0;
    z=M*r;
    p=z;
    rho_new=real(z'*r);
    for niter=1:n
        q=A*p;
        beta=real(p'*q);
        if beta/(q'*q)<tol2
            if beta<=0
                disp('Matrix A is not positive definite!');
                flag=0;
                break;
            end
            disp('Matrix A is ill conditioned');
        end
        rho=rho_new;
        alpha=rho/beta;
        dx=dx+alpha*p;
        r=r-alpha*q;
        norm_r2=r'*r;
        if sqrt(norm_r2)<tol
            flag=0;
            break;
        end
        z=M*r;
        rho_new=real(z'*r);
        rho_norm=rho_new/norm_r2;
        if rho_norm<tol2
            cur_iter=niter+waste_iter;
            if abs(rho_norm)>tol2
                fprintf('M is indefinite,      restarting. Iteration = %d \n',cur_iter);
            else
                fprintf('M is nearly singular, restarting. Iteration = %d \n',cur_iter);
            end
            gamma = delta*(tol2-rho_norm);
            waste_iter=waste_iter+niter;
            M=M + gamma*speye(n);
            x_0=x_0+dx;
            count=count+1;
            break;
        end
        p=z+(rho_new/rho)*p;
    end
end
x_0=x_0+dx;
toc;
iter=niter+waste_iter;
res=A*x_0-b;
err=norm(res);
disp(count);
end
