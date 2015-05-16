function b = b_kse (u)
n=8;
uu=zeros(1,24);
alpha=33;
uu(1:n) = u;
b=zeros(1,8);
for k=1:n
    mysum=0;
    for j=1:n
        if (k > j)
            kjsgn=1;
            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)));
        elseif (k < j)
            kjsgn=-1;
            mysum=mysum+j*u(j)*(uu(j+k)+kjsgn*uu(abs(k-j)));
        elseif (k == j)
            kjsgn=0;
            mysum=mysum+j*u(j)*uu(j+k);
        endif
    endfor
    b(k)=alpha*mysum/2;
endfor
endfunction
