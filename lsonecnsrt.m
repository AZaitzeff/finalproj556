function z=lsonecnsrt(D,y,L,z)
m=size(z,1);
z0=zeros(1,2*m);
z0(1:m)=z.*(z>0);
z0(m+1:2*m)=-z.*(z<0);
z = lsqlin([D -D],y,ones(1,2*m),L,[],[],zeros(1,2*m),[],z0);