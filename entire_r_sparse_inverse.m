function [M]=entire_r_sparse_inverse(A,n,lfil)
M=spalloc(n,n,lfil*n);
for j=1:n
    m=spalloc(n,1,lfil);
    r=sparse(j,1,1,n,1);
    for k=1:2*lfil
        [~, i]=max(abs(r));
        %del=r(i)/A(i,i); %this is modified for a scaled matrix
        del=r(i);
        m(i)=m(i)+del;
        if nnz(m)>=lfil
            break;
        end
        r=r-del*A(:,i);
    end
    M(:,j)=m;
end
end