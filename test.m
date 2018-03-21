rng(100)
factor=.5;
K=ceil(100*factor);
L=ceil(10*factor);
D=ceil(121*factor);
N=ceil(100*100*factor);
p=.01;

dt=randn(D,K);
for i=1:K
    dt(:,i)=dt(:,i)/(norm(dt(:,i))+.1);
end
zt=1.*(rand(L,N,K)<p);
y=zeros(N,L);
for i=1:L
   temp=zeros(N,1);
   for j=1:K
       dpad=[dt(:,j)' zeros(1,N-D)];
       Dm=toeplitz(dpad,[dpad(1) fliplr(dpad(2:end))]);
       temp=temp+Dm*reshape(zt(i,:,j),[N 1]);
   end
   y(:,i)=temp;
end

d0=randn(D,K);
for i=1:K
    d0(:,i)=d0(:,i)/(norm(d0(:,i))+.1);
end
z0=1.*(rand(L,N,K)<p);
%tic;[z,d,val,t]=skecting(y,D,K,.2,z0,d0,(p*K*N+10),(p*K*N+10));toc