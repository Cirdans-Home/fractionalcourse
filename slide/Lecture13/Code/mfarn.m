function [yout,err,flagsol,H,v]=mfarn(f,A,b,m,errt,rat);



%
%  errt    � la tolleranza per l'errore (espresso in log10)
%  rat     0 polinomiale 1 razionale
%
%  yout    approssimazione finale
%  err     vettore delle norme degli errori
%  flagsol segnalatore del raggiungimento della soluzione del sistema
%


tic;
flagsol=0;
[N,N]=size(A);
nr=norm(b);
v(:,1)=b/nr;

for k=1:m
    proj=zeros(N,1);
    if rat==0
        z=A*v(:,k);
    else
        z=A\v(:,k);
    end
    for l=1:k
        H(l,k)=v(:,l)'*z;
        proj=proj+H(l,k)*v(:,l);
    end
    vv=z-proj;
    H(k+1,k)=norm(vv);
    if H(k+1,k)~=0
        v(:,k+1)=vv/H(k+1,k);
    end
    tmp1=funm2(H(1:k,:),f);
    tmp2=nr*(tmp1(:,1));
    yout(:,k)=v(:,1:k)*tmp2;
    err(k)=abs(H(k+1,k))*abs(tmp2(k));
    if err(k)<errt
        flagsol=1;
        time = toc;
        fprintf('%d & %1.2e & %1.1e\n',k,err(k),time);
        %yout = yout(:,end);
        break;
    end
end
if flagsol == 0
    %yout = yout(:,end);
    time = toc;
    %fprintf('%d & %1.2e & %1.1e',k,err(k),time)
end
end

