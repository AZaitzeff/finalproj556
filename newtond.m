function x=newtond(A,b,t,x0,D)
[m,n]=size(A);
K=n/D;

ALPHA = 0.1; 
BETA = 0.5; 
MAXITERS = 1000; 
NTTOL = 5e-6; 
x=x0;
phi=zeros(1,K);
vals=zeros(1,K);
Hessb= t*(A'*A);
for iter = 1:MAXITERS
    grad = t*2*(A'*A*x-A'*b);
    Hess= Hessb;
    for i=1:K
        phi(i)=-log(1-sum(x((i-1)*D+1:i*D).^2));
        grad((i-1)*D+1:i*D)=grad((i-1)*D+1:i*D)+2*phi(i)*x((i-1)*D+1:i*D);
        xd=x((i-1)*D+1:i*D);
        Hess((i-1)*D+1:i*D,(i-1)*D+1:i*D)=Hess((i-1)*D+1:i*D,(i-1)*D+1:i*D)+4*(xd*xd')/phi(i)^2+2/phi(i)*eye(D);
    end
    v = -Hess\grad;
    lambda = grad'*v; 
    if abs(lambda)/2 < NTTOL, break; end;
    step = 1; 
    for i=1:K
        vals(i)=1-sum((x((i-1)*D+1:i*D)+step*v((i-1)*D+1:i*D)).^2);
    end
    
    while any(vals<=0) %backtrack until in domain
        step = BETA*step; 
        for i=1:K
            vals(i)=1-sum((x((i-1)*D+1:i*D)+step*v((i-1)*D+1:i*D)).^2);
        end
    end
    
    %while ( c'*(x+step*v)-sum(log(b-A*(x+step*v))) > val + ALPHA*step*lambda )  
    %    t = BETA*t; 
    %end
    x = x+step*v;
    
end