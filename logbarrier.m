function d=logbarrier(Z,y,d,mu,D,t0)
%Solves ||y-Zd||_2^2 with ||d_k||_2^2 <=1 where each d_k is dimension D. set mu =100 and t0=100 
%and if that does not work try mu =10 and t0=10
eps=.001;
t=t0;
MAXITER=1000;
K=size(d,1)/D;
for i=1:MAXITER
    d=newtond(Z,y,t,d,D);
    t=mu*t;
    if K/t<eps
        break
    end
end