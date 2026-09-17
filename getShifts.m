function [lambda2,mu2]=getShifts(l,d,relresBound,method)
% abs(l): dimension after restart，d: current dimension
% relresBound: 从端部向里排
global B B_bar
[~,~,~,C_tilde,S_tilde]=gsvd(B(1:d+1,1:d),B_bar(1:d,1:d),0); % gsvd升序
c=diag(C_tilde); s=diag(S_tilde);
if l>0
    c=flip(c); s=flip(s); % l>0: 从大到小；l<0: 从小到大
end

if method=="IRJBD"
    lambda2=c(abs(l)+1:d).^2;
elseif method=="IRRJBD"
    W=zeros(d,abs(l)); QA=zeros(d+1,abs(l)); QL=zeros(d,abs(l));
    for i=1:abs(l)
        tmp=[s(i)^2*B(1:d+1,1:d)'*B(1:d+1,1:d)-c(i)^2*B_bar(1:d,1:d)'*B_bar(1:d,1:d);...
            zeros(1,d-1),B(d+1,d+1)*B(d+1,d)];
        [~,~,Vtmp] = svd(tmp,'econ');
        W(:,i) = Vtmp(:,end);
        QA(:,i)=normalize(B(1:d+1,1:d)*W(:,i),"norm");
        QL(:,i)=normalize(B_bar(1:d,1:d)*W(:,i),"norm");
    end
    [Wfull,~]=qr(W);
    [QAfull,~]=qr(QA);
    [QLfull,~]=qr(QL);
    [~,~,~,C,S]=gsvd(QAfull(:,abs(l)+1:end)'*B(1:d+1,1:d)*Wfull(:,abs(l)+1:end),...
        QLfull(:,abs(l)+1:end)'*B_bar(1:d,1:d)*Wfull(:,abs(l)+1:end),0);
    lambda2=diag(C).^2./(diag(C).^2+diag(S).^2);
end

% 结合 adaptive shifting strategy，l含+adjust，relresBound不含
relgap=abs(c(abs(l))-relresBound(end)-sqrt(lambda2))/c(abs(l));
lambda2(relgap<1e-3)=double(l <= 0);
mu2=1-lambda2;
end