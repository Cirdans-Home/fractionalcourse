function v=dbprec(v)
%%DBPREC preconditioner for BVM system
%function v=dbprec(m,s,p,v,h_tilde,L,metodo,k)
%
% effettua il precondizionamento direttamente sui vettori v \tilde v
% con P=\breve A \otimes I_m-H*(B \otimes L)
%

global h_tilde hi L sqL m p s normL phi psi pol_flag OMEGA
global index1 index2 L1 U1 Lg Ug tmp onceLU Jc apprJ precond
%% AGGiUNTO DA FABIO PER PRECONDIZIONATORE AGGIUNTIVO
global FPR fracprec diagNdiv
%%

%% Operazioni aggiuntive per approccio flessibile:
if fracprec == 5 || fracprec == 7
    beta = bandwidth(FPR);
    omegaval=-1;
    OMEGAMAT = omegaval.^((0:m+beta-2)'/(m+beta-1));
    Lambdainv = zeros(m+beta-1,1);
end
%% Codice originario:
Im = speye(m);
sp=s+1; mp=m+1; ms=m*s; m2=m*2; mpp=m*p; flag=100;

% calcolo \tilde v=N'( (I_m \otimes F)Nv)
% N matrice di permutazione

v=v(index1); % permutaz.

for ipi=0:mpp-1
    i1=1+ipi*s; i2=(ipi+1)*s; % +++
    if precond==6
        v(i1:i2)=fft(conj(OMEGA.').*v(i1:i2),s);
    else
        v(i1:i2)=fft(v(i1:i2),s);
    end %if
end %for_ipi

v=v(index2); % antipermutaz.
j1=1; % j1 indice delle s1 fattorizzazioni LU necessarie (s1<=s)

% risolvo gli sp sistemi lineari di dimensione m con matrici G_j^i,
% i=1,...,p, j=1,...,s (lascia i primi s elementi inalterati)
for ipi=1:p
    %hi=h_tilde(ipi)/s;
    for j=1:s
        i1=1+m*(j-1)+ms*(ipi-1); i2=m*j+ms*(ipi-1); % +++
        v_tmp=v(i1:i2);
        if pol_flag
            flag=abs(tmp(j))*normL;
        end %if
        if ( flag<=1 && pol_flag ) % se usare approx polinomiale per (G_j^i)^-1
            v_tmp = (1/phi(j))*(speye(m)+L*tmp(j))*v(i1:i2); % appr. I ordine
        else
            if (apprJ==-1) % uso appross. tipo -1 per J (Aggiunto da Fabio: ovvero BCCB)
                %Jc_tmp=-hi*psi(j)*Jc; Jc_tmp=Jc_tmp+phi(j).*ones(m,1);
                %Jc_tmp=1./Jc_tmp;
                %Jc_tmp=1./(phi(j).*ones(m,1)-hi*psi(j)*Jc);
                %v_tmp=fft( Jc_tmp.*ifft(v_tmp,m),m);
                v_tmp=fft( ifft(v_tmp,m)./(phi(j).*ones(m,1)-hi*psi(j)*Jc),m);
                %fprintf('%e\n',norm(imag(v_tmp)));
                %                 cc = [L(1:ceil(m/2),1);L(1,m-ceil(m/2)+1:-1:2)'];
                %                 CL = gallery('circul',cc);
                %                 v_tmp= (phi(j)*speye(m)-hi*psi(j)*CL) \ v_tmp; % risolve tramite Gauss
            else
                if (onceLU)
                    i3=1+(j1-1)*m; i4=j1*m; j1=j1+1;
                    v_tmp = Ug(i3:i4,:)\(Lg(i3:i4,:)\v_tmp); % risolve tramite fatt. PLU
                elseif fracprec == 3 || fracprec == 4
                    v_tmp= (phi(j)*speye(m)-hi*psi(j)*FPR) \ v_tmp;
                    %Kj = (phi(j)*speye(m)-hi*psi(j)*FPR);
                    %v_tmp = tsolve(phi(j)*Im(:,1)-hi*psi(j)*FPR(:,1),phi(j)*Im(1,:)-hi*psi(j)*FPR(1,:),v_tmp);
                elseif fracprec == 5
                    Kj = (phi(j)*speye(m)-hi*psi(j)*FPR);
                    %c = [Kj(1:m,1);omegaval*Kj(1,beta:-1:2)'];
                    %Lambda = real(fft(full(OMEGAMAT.*c)));
                    %Lambdainv(Lambda>0)=1./Lambda(Lambda>0);
                    %Lambdainv(Lambda<=0)=0;
                    %[v_tmp,~] = gmres( Kj,v_tmp,20,1e-3,m, @(x) hankenagy( Lambdainv,x,m,beta,omegaval ));
                    %tic;
                    %[ preconvec ] = jacksonprec( full(Kj(1,:)),2 );
                    %[v_tmp] = gmres( Kj,v_tmp,20,1e-3,m,@(x) cprec(preconvec,x) );
                    %toc;
                    [v_tmp,~] = bicgstab(Kj,v_tmp,1e-3,20);
                    %                     [v_tmp,~] = gmres( Kj,v_tmp,20,1e-3,m);
                elseif fracprec == 6
                    % Precondizionatore Strang 2D per lo Jacobiano
                    SC = (phi(j)*speye(m)-hi*psi(j)*FPR);
                    c = reshape(SC(:, 1), sqrt(m) , []);     %// First column of each block
                    v_tmp = reshape(ifft2(fft2(reshape(v_tmp,sqrt(m),[]))./fft2(full(c))),m,1);
                elseif fracprec == 7
                    Kj = (phi(j)*speye(m)-hi*psi(j)*FPR);
                    [v_tmp,~] = gcrodr( Kj,v_tmp,20,5,[],1e-6);
                else % Precondizionatore con Jacobiano pieno
                    %% Risolve tramite Gauss
                    v_tmp= (phi(j)*speye(m)-hi*psi(j)*L) \ v_tmp; 
                    %% Risolve tramite Toeplitz Solver:
                    %v_tmp=
                    %tsolve(phi(j)*Im(:,1)-hi*psi(j)*L(:,1),phi(j)*Im(1,:)-hi*psi(j)*L(1,:),v_tmp);
                    
                end %if onceLU
            end %if apprJ
        end %if flag<=1
        v(i1:i2)=v_tmp;
    end %j
end %ipi

% calcolo il vettore finale 'precondizionato'
% r=N'(I_m \otimes F^*)N \tilde v)
% N matrice di permutazione

v=v(index1); % permutaz.

for ipi=0:mpp-1
    i1=1+ipi*s; i2=(ipi+1)*s; % +++
    v(i1:i2)=ifft(v(i1:i2),s);
    if precond==6
        v(i1:i2)=(OMEGA.').*v(i1:i2);
    end %if
end %for_ipi

v=v(index2); % antipermutaz. finale

v=real(v); % le parti immaginarie sono dovute ad errori di arrotondam.

% end dbprec.m

