function imRestartJBD(l,k,lambda2,mu2)
% Implicitly restart k-step JBD process to l-step. The shifts satisfies lambda2+mu2=1
global U U_hat V_prime B B_bar;
g=[zeros(1,k),1]; %e_{k+1}'*G
for i=1:k-l
    [G,~]=givrot([B(1,1)^2-lambda2(i),B(1,1)*B(2,1)]); % corresponding to G_1^T in the paper
    B(1:2,1:2)=G*B(1:2,1:2);
    if B(1,1)<0
        G=-G;
        B(1:2,1:2)=-B(1:2,1:2);
    end
    U(:,1:2)=U(:,1:2)*G';
    g(1:2)=g(1:2)*G';
    
    [J,r]=givrot(B(1,1:2)); % J_1^T
    B(1,1:2)=[r,0];
    B(2:3,1:2)=B(2:3,1:2)*J';
    B_bar(:,1:2)=B_bar(:,1:2)*J';
    V_prime(:,1:2)=V_prime(:,1:2)*J';
    
    [G_bar,~]=givrot(B_bar(1:2,1)); % \bar G_1^T
    B_bar(1:2,:)=G_bar*B_bar(1:2,:);
    U_hat(:,1:2)=U_hat(:,1:2)*G_bar';
    for j=2:k-1 % J_j and G_j  
        [G,r]=givrot(B(j:j+1,j-1)); % G_j^T
        B(j:j+1,j-1)=[r;0];
        B(j:j+1,j:j+1)=G*B(j:j+1,j:j+1);
        U(:,j:j+1)=U(:,j:j+1)*G';
        g(j:j+1)=g(j:j+1)*G';
        
        [J,r]=givrot(B(j,j:j+1)); % J_j^T
        B(j,j:j+1)=[r,0];
        B(j+1:j+2,j:j+1)=B(j+1:j+2,j:j+1)*J';
        B_bar(:,j:j+1)=B_bar(:,j:j+1)*J';
        V_prime(:,j:j+1)=V_prime(:,j:j+1)*J';
        
        [G_bar,~]=givrot(B_bar(j:j+1,j)); % \bar G_j^T
        B_bar(j:j+1,1:k)=G_bar*B_bar(j:j+1,1:k);
        U_hat(:,j:j+1)=U_hat(:,j:j+1)*G_bar';
    end
    [G,r]=givrot(B(k:k+1,k-1)); % G_k^T
    B(k:k+1,k-1)=[r;0];
    B(k:k+1,k)=G*B(k:k+1,k);
    U(:,k:k+1)=U(:,k:k+1)*G';
    g(k:k+1)=g(k:k+1)*G';
end
[V_prime(:,l+1),~,B(l+1,l+1)]=normalize(B(k+1,k+1)*g(l+1)*V_prime(:,k+1)+B(l+1,l+1)*V_prime(:,l+1),"norm");
B_bar(l,l+1)=-B(l+1,l+1)*B(l+1,l)/B_bar(l,l);
B_bar=[diag(diag(B_bar))+diag(diag(B_bar(:,1:k),1),1),[zeros(k-1,1);B_bar(k,k+1)]];
end

function [G,r]=givrot(a) % Numerical Methods for Least Squares Problems p54 Algorithm2.3.2
x=a(1); y=a(2); % G*[x;y]=[r;0]
if y==0
    c=1.0; s=0.0; r=x;
elseif abs(y)>abs(x)
        t=x/y; tt=sqrt(1+t^2);
        s=1/tt; c=t*s; r=tt*y;
else
    t=y/x; tt=sqrt(1+t^2);
    c=1/tt; s=t*c; r=tt*x;
end
G=[c,s;-s,c];
end