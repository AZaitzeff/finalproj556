function z=lassocd(D,y,lam,z0)
    S=@(z,gam)(sign(z)*max(abs(z)-gam,0));
    [N,p]=size(D);
    tol=1e-4;
    MAXITER=1000;
    z=z0;
    yt=D*z;
    %inner=D'*y;
    %r=y-yt;
    %Z = sparse(p,p);
    %index=1:p;
    %val=zeros(1,p);
    val=vecnorm(D).^2;
    %for i=1:p
    %    val(i)=D(:,i)'*D(:,i);
    %end
    for i=1:p
        %if z(i)~=0
        %    yt=yt-D(:,i)*z(i);
        %end
        z(i)=S(D(:,i)'*(y-yt)/val(i),lam/val(i));
        if z(i)~=0
            yt=yt+D(:,i)*z(i);
        end
    end
    index=find(z)';
    for t=1:MAXITER
        for i=index
            yt=yt-D(:,i)*z(i);
            z(i)=S(D(:,i)'*(y-yt)/val(i),lam/val(i));
            yt=yt+D(:,i)*z(i);
            %z(i)=S((D(:,i)'*r+val(i)*z(i))/val(i),lam/val(i));
            %if abs(z(i))>tol
            %    r=y-D*z;
            %end
            %temp=S((inner(i)-Z(i,:)*z+val(i)*z(i))/val(i),lam/val(i));
            %if abs(temp)>tol && abs(z(i))<tol
            %    Z(:,i)=D(:,i)'*D;
            %end
            %z(i)=temp;
        end
        
        sum((z(index)-z0(index)).^2);
        if sum((z-z0).^2)<tol*N
            break;
        end
        index=find(z)';
        z0=z;
    end
end

