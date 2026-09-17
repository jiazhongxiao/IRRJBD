function testJBD(d)
% compute and print errors of current JBD process (d-step)
global U U_hat V_prime B B_bar;
m=size(U,1); p=size(U_hat,1);
tmpV_prime=V_prime(:,1:d); tmpU=U(:,1:d+1); tmpU_hat=U_hat(:,1:d);
tmpB=B(1:d+1,1:d); tmpB_bar=B_bar(1:d,1:d);

fprintf("Orthogonality of V_prime, U, U_hat: %e, %e, %e\n",norm(tmpV_prime'*tmpV_prime-eye(d)),...,
   norm(tmpU'*tmpU-eye(d+1)),norm(tmpU_hat'*tmpU_hat-eye(d)));

fprintf("inv of B(1:end-1,:) and B_bar: %e, %e\n",norm(inv(tmpB(1:d,:))),norm(inv(tmpB_bar)));
fprintf("minim of alpha, beta: %e\n",min([min(abs(diag(tmpB))),min(abs(diag(tmpB,-1)))]));
fprintf("minim of alpha_hat, beta_hat: %e\n",min([min(abs(diag(tmpB_bar))),min(abs(diag(tmpB_bar,1)))]));
fprintf("error of BTB+BbarTBbar=I: %e\n",norm(tmpB'*tmpB+tmpB_bar'*tmpB_bar-eye(d)));
end