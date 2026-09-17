function relresBound=getrelresBound(target,d,method)
% relresBound是从端部到内部的
global B B_bar
[QA_tilde,QL_tilde,~,C_tilde,S_tilde]=gsvd(B(1:d+1,1:d),B_bar(1:d,1:d),0);
c=diag(C_tilde); s=diag(S_tilde);
if target>0 % target>0则大的在前；target<0则小的在前
    QA_tilde=flip(QA_tilde,2); QL_tilde=flip(QL_tilde,2); c=flip(c); s=flip(s);
end

if method=="IRJBD"
    relresBound=abs(B(d+1,d+1)*QA_tilde(d,1:abs(target))*diag(s(1:abs(target))) ...
        -B_bar(d,d+1)*QL_tilde(d,1:abs(target))*diag(c(1:abs(target))));
elseif method=="IRRJBD"
    relresBound=zeros(abs(target),1);
    for i=1:abs(target)
        tmp=[s(i)^2*B(1:d+1,1:d)'*B(1:d+1,1:d)-c(i)^2*B_bar(1:d,1:d)'*B_bar(1:d,1:d);...
            zeros(1,d-1),B(d+1,d+1)*B(d+1,d)];
        [~,~,Vtmp] = svd(tmp,'econ');
        w = Vtmp(:,end);
        c_hat=norm(B(1:d+1,1:d)*w);
        s_hat=norm(B_bar(1:d,1:d)*w);
        relresBound(i)=sqrt(norm(s_hat^2*B(1:d+1,1:d)'*B(1:d+1,1:d)*w ...
            -c_hat^2*B_bar(1:d,1:d)'*B_bar(1:d,1:d)*w)^2 ...
            +B(d+1,d+1)^2*B(d+1,d)^2*w(end)^2)/(c_hat*s_hat);
    end
end
end