addpath('IRRJBDtools');
global A Aname L Lname U U_hat V_prime B B_bar;
files = dir("../IRRJBDmatrices/*.mat");
Anames = erase(string({files.name}), ".mat"); 
Atrans=[0]; % whether to transpose A
Lnames=["Lflat"];
matrixInfoPath="";
resultPath="../results.csv";
targets=[-5,-10];
ks=[25,50];
adjust=3;
tol=1e-14;
maxit=5000;
reorth=1; % 0: no reorthogonalization，1: full reorthogonalization 2: reorthogonalization only on V_prime and U, not on U_hat
lsqrtol=10*eps;
lsqrmaxitcoef=10; % lsqrmaxit=lsqrmaxitcoef*n
argCombs=allcomb(Anames,Atrans,Lnames,targets,ks);
argCombsCell=num2cell(argCombs);

for argIdx=1:size(argCombsCell,1)
    [Aname,Atrans,Lname,target,k]=argCombsCell{argIdx,:};
    Atrans=str2num(Atrans); target=str2num(target); k=str2num(k);
    [m,n,p]=loadmatrix(sprintf("../IRRJBDmatrices/%s.mat",Aname),Atrans,Lname,matrixInfoPath);
    if Atrans==1
        Aname=Aname+"^T";
    end
    lsqrmaxit=round(lsqrmaxitcoef*n);
    fprintf("Aname: %s, Lname: %s, target: %d, k: %d\n",Aname,Lname,target,k);
    
    rng(2024); % random seed
    u1=normalize(randn(m,1),"norm");
    U=zeros(m,k+1); U(:,1)=u1; U_hat=zeros(p,k); V_prime=zeros(m+p,k+1);
    B=zeros(k+1,k+1); B_bar=zeros(k,k+1);
    t0=cputime;
    [d,Res_b1,relresVec,flag1]=IRJBD(target,k,adjust,reorth,tol,maxit,lsqrtol,lsqrmaxit,"IRJBD");
    time1=cputime-t0;
    relres=evaluate(target,d,lsqrtol,lsqrmaxit);
    iter1=size(relresVec,1); Res1=relres; 
    invBB_bar1=norm(inv(B(1:d,1:d)))*norm(inv(B_bar(1:d,1:d)));
    save(sprintf("../relresVec/%s,%s,IR,%d,%d.mat",Aname,Lname,target,k),"relresVec");
    
    U=zeros(m,k+1); U(:,1)=u1; U_hat=zeros(p,k); V_prime=zeros(m+p,k+1);
    B=zeros(k+1,k+1); B_bar=zeros(k,k+1);
    t0=cputime;
    [d,Res_b2,relresVec,flag2]=IRJBD(target,k,adjust,reorth,tol,maxit,lsqrtol,lsqrmaxit,"IRRJBD");
    time2=cputime-t0;
    relres=evaluate(target,d,lsqrtol,lsqrmaxit);
    iter2=size(relresVec,1); Res2=relres; 
    invBB_bar2=norm(inv(B(1:d,1:d)))*norm(inv(B_bar(1:d,1:d)));
    save(sprintf("../relresVec/%s,%s,IRR,%d,%d.mat",Aname,Lname,target,k),"relresVec");

    if resultPath~=""
        variablenames={'A','L','target','k','iter1','iter2','SI(\%)','time1','time2','ST(\%)',...,
            'Res_b1','Res_b2','Res1','Res2','invB*invB_bar1','invB*invB_bar2','flag1','flag2'};
        T=table(Aname,Lname,target,k,iter1,iter2,(iter1-iter2)/iter1*100,time1,time2,(time1-time2)/time1*100,...,
            Res_b1,Res_b2,Res1,Res2,invBB_bar1,invBB_bar2,flag1,flag2,'VariableNames',variablenames);
        writetable(T,resultPath,'WriteMode','append');
    end
end
rmpath('IRRJBDtools');