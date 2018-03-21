function [z,d,val,t]=skecting(y,D,K,alpha,z0,d0,m,s)
%TO DO right now it only works for vectors. Needs to implement skecting
[N,L]=size(y);
%z=zeros(L,N,K);
%d=zeros(D,K);
z=z0;
d=d0;
val0=inf;
%dflat=reshape(d,[D*K 1]);
%zflat=reshape(z(1,:,:),[N*K 1]);
tol=1e-5;
MAXITER=1000;
for t=1:MAXITER
    Ad=zeros(N,N*K);
    for i=1:K
        dpad=[d(:,i)' zeros(1,N-D)];
        Ad(:,(i-1)*N+1:i*N)=toeplitz(dpad,[dpad(1) fliplr(dpad(2:end))]);
    end
    P = sjlt(m, N, s);
    for i=1:L
        zflat=reshape(z(i,:,:),[N*K 1]);
        zflat=lassocd(Ad,y(:,i),alpha,zflat);
        %zflat=lsonecnsrt(Ad,y(:,i),alpha,zflat);%solves ||y-Dz||^2_2 with constraint on z.
        z(i,:,:)=reshape(zflat, [N K]);
    end
    Az=zeros(N*L,D*K);%convolution Matrix z
    for j=1:L
        for i=1:K
            ztemp=reshape(z(j,:,i),[1 N]);
            toeptemp=toeplitz(ztemp,[ztemp(1) fliplr(ztemp(2:end))]);
            Az((j-1)*N+1:j*N,(i-1)*D+1:i*D)=toeptemp(:,1:D);
        end
    end
    dflat=reshape(d,[D*K 1]);
    yflat=reshape(y,[N*L 1]);
    m2=ceil(1.05*D*K);
    %s2=ceil(log(N*L)^2);
    P = sjlt(m2, N*L, D);
    dflat=logbarrier(Az,yflat,dflat,10,D);%solves ||y-Zd||^2_2 with constraint on d. uses barrier method that I implemented
    d=reshape(dflat,[D K]);
    val=sum((yflat-Az*dflat).^2)+alpha*sum(abs(z(:)))
    chng=sum((d(:)-d0(:)).^2)+sum((z(:)-z0(:)).^2)
    if chng<tol*N || val>val0
        break;
    end
    d0=d;
    z0=z;
    val0=val;
end

