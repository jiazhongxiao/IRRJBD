function [m,n,p]=loadmatrix(Apath,Atrans,Lname,csvpath)
% load matrices A and L and save the informations to csvpath (if provided)
% Lname: L1: diag(3,1,1); L2: diag(1,1) with p=n+1
global A L
S=load(Apath);
A=S.Problem.A;
if Atrans==1
    A=A';
end
m=size(A,1); n=size(A,2);

if Lname=="Ltall"
    i=[1:n,2:n+1];
    j=[1:n,1:n];
    v=[2*ones(1,n),ones(1,n)];
    p=n+1;
elseif Lname=="Lflat"
    i=[1:n-1,1:n-1];
    j=[1:n-1,2:n];
    v=[2*ones(1,n-1),ones(1,n-1)];
    p=n-1;
end
L=sparse(i,j,v,p,n);

if csvpath~=""
    [~, Aname] = fileparts(Apath);% 去掉路径和扩展名
    if Atrans==1
        Aname=Aname+"^T";
    end
    kappa1=0;
    sigmamax=svds([A;L],1,'largest','SubspaceDimension',1000);
    sigmamin=svds([A;L],1,'smallest','SubspaceDimension',1000);
    if ~isnan(sigmamin) && sigmamin>0
        kappa2=sigmamax/sigmamin;
    else
        kappa2=Inf;
    end
    variablenames={'A','L','m','n','p','nnz','kappa(A)','kappa([A;L])'};
    T=table(Aname,Lname,m,n,p,nnz(A)+nnz(L),kappa1,kappa2,'VariableNames',variablenames);
    writetable(T,csvpath,'WriteMode','append');
end
end