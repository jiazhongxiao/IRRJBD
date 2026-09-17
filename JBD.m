function [d,relresBound,flag,relres]=JBD(target,d,k,reorth,tol,lsqrtol,lsqrmaxit,method)
% Input: l: number of targets, l>0 stands for largest, while l<0 stands for smallest
%        d: current step of the JBD process，k: target step
% Output: Expand the JBD process to k-step.
%         d: step of JBD. If it converges early, it is possible that d<k
%         flag=0: converged，flag=1: not converged
global A L U U_hat V_prime B B_bar;
m=size(A,1); p=size(L,1);
flag=1;
if d==0
    [x,lsqrflag]=lsqr([A;L],[U(:,1);zeros(p,1)],lsqrtol,lsqrmaxit);
    if lsqrflag~=0
        fprintf("flag of lsqr=%d in %dth column\n",lsqrflag,1);
    end
    [V_prime(:,1),~,B(1,1)]=normalize([A;L]*x,"norm");
end

for i=d+1:k
    %扩张到i维
    if i==1
        [U_hat(:,1),~,B_bar(1,1)]=normalize(V_prime(m+1:m+p,1),"norm");
    else
        tmp=(-1)^(i-1)*V_prime(m+1:m+p,i)-(-1)^(i-1)*B_bar(i-1,i)*U_hat(:,i-1);
        if reorth==1
            tmp=tmp-U_hat(:,1:i-1)*(U_hat(:,1:i-1)'*tmp);
        end
        [U_hat(:,i),~,B_bar(i,i)]=normalize(tmp,"norm");
    end
    B_bar(i,i)=(-1)^(i-1)*B_bar(i,i);
    tmp=V_prime(1:m,i)-B(i,i)*U(:,i);
    if reorth~=0
        tmp=tmp-U(:,1:i)*(U(:,1:i)'*tmp);
    end
    [U(:,i+1),~,B(i+1,i)]=normalize(tmp,"norm");
    [x,lsqrflag]=lsqr([A;L],[U(:,i+1);zeros(p,1)],lsqrtol,lsqrmaxit);
    if lsqrflag~=0
        fprintf("flag of lsqr=%d in %dth column\n",lsqrflag,i+1);
    end
    tmp=[A;L]*x-B(i+1,i)*V_prime(:,i);
    if reorth~=0
        tmp=tmp-V_prime(:,1:i)*(V_prime(:,1:i)'*tmp);
    end
    [V_prime(:,i+1),~,B(i+1,i+1)]=normalize(tmp,"norm");
    B_bar(i,i+1)=-B(i+1,i+1)*B(i+1,i)/B_bar(i,i);
    
    if i>=abs(target)
        relresBound=getrelresBound(target,i,method);
        if max(relresBound)<tol
            flag=0;
            break
        end
    end
end
d=i;
relres=0;
%relres=evaluate(target,d,10*eps,lsqrmaxit);
end