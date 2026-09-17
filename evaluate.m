function relres=evaluate(target,d,lsqrtol,lsqrmaxit)
global A L U U_hat V_prime B B_bar;
Rbound=sqrt(norm(A,1)*norm(A,inf)+norm(L,1)*norm(L,inf));

[~,~,~,C,S]=gsvd(B(1:d+1,1:d),B_bar(1:d,1:d),0);
% C=svd(B(1:d+1,1:d),0); C=diag(C(end:-1:1));
% S=svd(B_bar(1:d,1:d)); S=diag(S);
X=zeros(size(A,2),abs(target)); Y=zeros(size(A,1),abs(target)); Z=zeros(size(L,1),abs(target));
if target>0
    C=C(d:-1:d-target+1,d:-1:d-target+1); CA=C;
    S=S(d:-1:d-target+1,d:-1:d-target+1); SL=S;
else
    C=C(1:-target,1:-target); CA=C;
    S=S(1:-target,1:-target); SL=S;
end
for i=1:abs(target)
    tmp=[S(i,i)^2*B(1:d+1,1:d)'*B(1:d+1,1:d)-C(i,i)^2*B_bar(1:d,1:d)'*B_bar(1:d,1:d);...
        zeros(1,d-1),B(d+1,d+1)*B(d+1,d)];
    %[~,~,w]=svds(tmp,1,"smallest");
    [~,~,Vtmp] = svd(tmp,'econ');
    w = Vtmp(:,end); %250923：svds会警告矩阵接近奇异，改成svd试试
    [X(:,i),~]=lsqr([A;L],V_prime(:,1:d)*w,lsqrtol,lsqrmaxit);
    [qA,~,CA(i,i)]=normalize(B(1:d+1,1:d)*w,"norm");
    Y(:,i)=U(:,1:d+1)*qA;
    [qL,~,SL(i,i)]=normalize(B_bar(1:d,1:d)*w,"norm");
    Z(:,i)=U_hat(:,1:d)*qL;
end

relres=0;
for i=1:abs(target)
    r1=norm(A*X(:,i)-CA(i,i)*Y(:,i));
    r2=norm(L*X(:,i)-SL(i,i)*Z(:,i));
    r3=norm(SL(i,i)*A'*Y(:,i)-CA(i,i)*L'*Z(:,i));
    relres=max(relres,sqrt(r1^2+r2^2+r3^2)/Rbound);
end
end