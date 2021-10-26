function [x,i,flag] = spcgls(A,b,tol,maxit,M,x0)
[m,n]=size(A);
if (nargin<3) tol=1e-8; end
if (nargin<4) maxit=min([m,n]); end
if (nargin<5) M=speye(n); end
if (nargin<6) x0=zeros(n,1); end
tol2=1e-2;
delta=10;
count=0;
At=A';
flag=2;
waste_iter=0;
while(flag==2)
    r=b-A*x0;
    x=spalloc(n,1,n);
    t=At*r;
    norm0=norm(t);
    tol=tol*norm0;
    w=M*t;
    gamman=t'*w;
    u=w;
    for i=1:maxit
        gamma=gamman;
        q=A*u;
        q2=q'*q;
        if (q2<=0)
            flag=-1; %indefinite system
            break;
        end
        alpha=gamma/q2;
        x=x+alpha*u;
        r=r-alpha*q;
        t=At*r;
        normt=norm(t);
        if (normt<tol)
            flag=0; %solved
            break
        end
        w=M*t;
        gamman=t'*w;
        if (gamman/normt^2)<tol2
            waste_iter=waste_iter+i;
            if abs(gamman)>tol2
                fprintf('M is indefinite,      restarting. Iteration = %d \n',waste_iter);
            else
                fprintf('M is nearly singular, restarting. Iteration = %d \n',waste_iter);
            end
            eta = delta*(tol2-gamman);
            M=M + eta*speye(n);
            x0=x0+x;
            count=count+1;
            break;
        end
        beta=gamman/gamma;
        u=w+beta*u;
        if (i==maxit) flag=1; end %no convergence
    end
end
i=i+waste_iter;
x=x0+x;
end