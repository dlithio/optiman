u = [ -2.05706390e-29,   5.72735089e+00,   1.10957111e-28,-1.96497363e+00,  -3.14574680e-29,   2.74222890e-01,1.11849717e-29,  -3.23815403e-02];

f_kse(u);

realu = fsolve (@f_kse,u);

jac = zeros(8,8);
h = 1e-6;
for j=1:8
forwardu = realu;
forwardu(j) = forwardu(j) + h;
backwardu = realu;
backwardu(j) = backwardu(j) - h;
forwardf = f_kse(forwardu);
backwardf = f_kse(backwardu);
jac(:,j) = (forwardf - backwardf)/(2*h);
endfor

[V, lambda] = eig (jac)
