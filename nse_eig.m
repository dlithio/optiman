jac = textread("jac_output_59");
jac = reshape(jac,112,112);
[V, lambda] = eig (jac);
idx = find(real(diag(lambda)) > 0);
eigvec1 = V(:,idx(1));
eigvec2 = V(:,idx(2));
save -ascii "eigvec1" eigvec1
save -ascii "eigvec2" eigvec2
