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

AA=A'*A;
Ab=A'*b;
Hessb= t*(AA);
%gradb = t*2*(A'*A*x-A'*b);
for iter = 1:MAXITERS
    xtemp=reshape(x,[D,K]);
    vals=sum(xtemp.^2,1);
    val=t*x'*(AA*x-2*A'*b)-sum(log(1-vals));
    grad = t*2*(AA*x-Ab);
    Hess= Hessb;
    phi=1-vals;
    for i=1:K
        %phi(i)=(1-sum(x((i-1)*D+1:i*D).^2));
        grad((i-1)*D+1:i*D)=grad((i-1)*D+1:i*D)+2/phi(i)*x((i-1)*D+1:i*D);
        xd=x((i-1)*D+1:i*D);
        Hess((i-1)*D+1:i*D,(i-1)*D+1:i*D)=Hess((i-1)*D+1:i*D,(i-1)*D+1:i*D)+4*(xd*xd')/phi(i)^2+2/phi(i)*eye(D);
    end
    %m2=ceil(n/4);
    %P = sjlt(m2, n, m2);
    v = -(Hess)\grad;
    %v = -grad./diag(Hess);
    lambda = grad'*v;
    if abs(lambda)/2 < NTTOL, break; end;
    step = 1;
    xv = x+step*v;
    xtemp=reshape(xv,[D,K]);
    vals=sum(xtemp.^2,1);
    
    while any(vals>1) %backtrack until in domain
        step = BETA*step; 
        xv = x+step*v;
        xtemp=reshape(xv,[D,K]);
        vals=sum(xtemp.^2,1);
    end
    
    while ( t*(xv'*AA*xv-2*xv'*Ab)-sum(log(1-vals)) > val + ALPHA*step*lambda )  
        step = BETA*step; 
        xv = x+step*v;
        xtemp=reshape(xv,[D,K]);
        vals=sum(xtemp.^2,1);
    end
    x = xv;
    
end