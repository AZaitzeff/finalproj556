rng(100)

K=2;
L=10;
D=4;
N=100;
p=.2;

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
[z,d,val,t]=skecting(y,D,K,p*K*N+10,z0,d0);