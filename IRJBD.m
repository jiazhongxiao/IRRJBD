function [d,relresBound,relresVec,flag]=IRJBD(target,k,adjust,reorth,tol,maxit,lsqrtol,lsqrmaxit,method)
global Aname Lname;
d=0; l=abs(target);
relresVec=zeros(maxit,1);
try
    for iter=1:maxit
        [d,relresBound,flag,relres]=JBD(target,d,k,reorth,tol,lsqrtol,lsqrmaxit,method);
        relresVec(iter)=relres;
        if flag==0
            break
        end
        [lambda2,mu2]=getShifts(sign(target)*(l+adjust),d,relresBound,method);
        imRestartJBD(l+adjust,k,lambda2,mu2);
        d=l+adjust;
        if iter==1 || mod(iter,10)==0
            testJBD(d);
            fprintf("%s, %s, %s, target=%d, k=%d, iter=%d, relresBound=%e\n",Aname,Lname,method,target,k,iter,max(relresBound));
        end
    end
catch
    flag=2;
end
relresVec=relresVec(1:iter);
relresBound=max(relresBound);
end