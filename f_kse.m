function f = f_kse (u)
n=8;
alpha=33;
b = b_kse (u);
f = zeros(1,8);
for j=1:n
    f(j)=(-4*j .^ 4)*u(j)+alpha*(j .^ 2)*u(j) - b(j);
endfor
endfunction
