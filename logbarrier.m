function d=logbarrier(Z,y,d,mu,D)
eps=.001;
t=1;
MAXITER=1000;
K=size(d,1)/D;
for i=1:MAXITER
    d=newtond(Z,y,t,d,D);
    t=mu*t;
    if K/t<eps
        break
    end
end