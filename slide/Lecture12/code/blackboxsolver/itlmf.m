% itlmf.m -  versione 1.5.2 - 21/12/2000
%
% Integratore per ODE, PDE in forma di
% IVP (lineari qui) con metodi multistep
% con condizioni iniziali (IVM) o iniz./bordo (BVM)
%
% Sono considerati metodi iterativi  precondizionati dalla
% nuova classe di matrici definite da D. Bertaccini
%
% codice by D. Bertaccini
%
% Universita' di Roma "La Sapienza"
% Dipartimento di Matematica "G. Castelnuovo"
% Piazzale A. Moro 2  00185 Roma, Italy
% tel. 06-4991-3253, fax: 06-447.01.007/+39-06-4991-3253
% home: 055-75.46.41, 347-84.68.559
% e-mail: bertaccini@mat.uniroma1.it
% http://www.mat.uniroma1.it/~bertaccini
%
% richiede: mab, rosi_gam, rosi_bdf, rosi_etr, rosi_et2, rosi_tom
% coef, vsolve, cvsolve, dbprecb, dbprec, bi_cg.m, bicgsta2
% gmres, ugmres, bicgspsx, bicgspdx, bicg2psx
%
% Modificato da Fabio Durastante per adattamenti alle equazioni frazionarie
% 21/12/2016
% Modificato da Fabio per girare su una suite di tanti test per un medesimo
% precondizionatore

clear
clear global
clc
%close all
%

global precond pol_flag fracprec apprJ diagNdiv

pb=input('Problema[1-20] (Articolo Originario) [21-25] (Frazionarie):');
if pb == 23 || pb == 24
    fracval = input('Esponente frazionario: ');
end
disp('0=no precond. 2=T.Chan 3=Strang 4=P-circulant 5=modified-Strang 6=skew 7=Omega 8=Strang-Jm')
disp('9=Skew-Jm 10=LowMem-Jm 11=LowMem-Jm Fixed Number of Diagonals 12=Strang-Jm-2D');
precond=input('Scelta [0-10] :');
if precond == 8
    apprJ = -1;
    %     disp('2=T.Chan 3=Strang 4=P-circulant 5=modified-Strang 6=skew 7=Omega')
    %     precond=input('Scelta [2-7] :');
    precond = 3; % Selezione Strang per A,B,Jm
    fracprec = 1;
elseif precond == 9
    apprJ = -1;
    %     disp('2=T.Chan 3=Strang 4=P-circulant 5=modified-Strang 6=skew 7=Omega')
    %     precond=input('Scelta [2-7] :');
    precond = 6; % Selezione Skew-Circolant per A,B,Jm
    fracprec = 2;
elseif precond == 10
    diagNdiv = input('Frazione delle diagonali da preservare: ');
    apprJ = 0;
    precond = 6; % Selezione Skew-Circolant per A,B e g_k(J_m)
    fracprec = 3;
elseif precond == 11
    diagNdiv = input('Numero delle diagonali da preservare: ');
    apprJ = 0;
    precond = 6; % Selezione Skew-Circolant per A,B e g_k(J_m)
    fracprec = 4;
elseif precond == 12
    apprJ = 0;
    precond = 3; % Selezione Strang per A,B,Jm
    fracprec = 6;
else
    apprJ = 0;
    fracprec = 0;
end
metodo=input('-1=MIDPOINT 0=GBDF 1=GAMs 2=ETR2 3=TOMs 4=ETR2/1:');
passi=input('k (numero passi LMF) =');

fprintf('Solver Lineare:\n');
fprintf('0)BiCGSTAB(2)\n1)GMRES(20)\n2) BiCGstab\n3)FGMRES+Variable DBPREC\n10) Preconditioned IDR(4)\n11) FGMRES+GCRODR\n')
mspl = input('Solver = ');
if mspl==3
    fracprec=5;
elseif mspl==11
    fracprec=7;
end

%for svalue = [64,128,256,512,1024,2048,4096,8192]/2%[25,49,97,193,385,769,1537]%[24,48,96]+1%[33,65,129,257]%
    for  Nvalue = [64,128,256,512,1024,2048,4096,8192]%[16,32,64,96,128]%
        svalue = Nvalue/2;
        %fprintf('\n------------------------------------------------------\n');
        %fprintf('N = %d, s = %d\n',Nvalue,svalue);
        %% SETTING OPZIONI

        global normL phi psi coe index1 index2 h_tilde hi residuo_vero
        global L sqL m p s Lsparsa onceLU nzL
        global OMEGA
        global FPR Jc % Flag e variabile per i precondizionatori per
        % Jacobiani frazionari
        global EV % Variable per il BCCB

        IVPfun=0; hi=-1; iter=0; % mspl=-1;
        non_prec_plot=0; % se ==0 non calcola eigensystem di M
        icn_yes=0; % calcola gli icn e sum(1/icn(1:m*s*p));
        % da usare con le molle se XX_p malcond.

        %plot_spectra=input('Plotto lo spettro di P^{-1}*M ? [0/1]:');
        plot_spectra=0;     % Calcola condizionamento di P, M e se eig_computation=1
        eig_computation=0;  % anche lo spettro degli autovalori
        flops=0;

        %pb=input('Problema[1-20] (Articolo Originario) [21-24] (Frazionarie):');
        %pb = 23;

        %precond=input('P:2=opt 3=Str 4=Dan 5=mStr 6=skew 7=Omega 8=Felix:');
        %         disp('Precondizionatore:')
        %         disp('0=no precond. 2=T.Chan 3=Strang 4=P-circulant 5=modified-Strang 6=skew 7=Omega 8=Strang-BCCB 9=Skew-LowMem')
        %         precond=input('Scelta [0-7] :');
        %        precond = 8;

        if precond~=0
            pol_flag=1; % pol_flag==1 permette approx. polinomiale
        else
            pol_flag=0; % Strang ha \breve A, \breve B singolari
        end %if
        pol_flag=0;

        %         metodo=input('-1=MIDPOINT 0=GBDF 1=GAMs 2=ETR2 3=TOMs 4=ETR2/1:');
        %         passi=input('k (numero passi LMF) =');
        %         metodo = 0;
        %         passi = 2;

        restart_it=20;  % GMRES(restart_it)

        if pb<0
            error('Problema sconosciuto!\n');
        end

        % fprintf('\n1) GMRES(%i)\n',restart_it);
        % disp('2) BiCGstab')
        % disp('3) BiCG')
        % disp('4) BiCGstab(2)')
        % disp('5) CG-normal eq')
        % disp('8) Circulant')
        % disp('9) sparse LU')
        % if mspl>=0, else, mspl=input(':'); end %if

        %mspl=1;

        if (mspl~=9)&&(mspl~=0)&&(mspl~=3)&&(mspl~=4)&&(mspl~=6)
            %onceLU: 0=fattorizzo G_j x iterazione, 1=una tantum+mem, -1=approssimo G_j con matrice circ.
            %onceLU=input('onceLU=[0/1]');
            onceLU=0;
        end %if mspl

        %tol=input('Tolleranza su ||r||:');
        tol=1e-8; % tolleranza su residuo relativo del metodo iterativo
        residuo_vero=1; % 1=tolleranza sul residuo vero r=b-M*x

        %s=input('s (dim. matrici A e B, consiglio potenza di 2)='); %
        s = svalue;

        if (0) % hi costante deciso a priori
            if (pb>9)
                if (s<20)
                    hi=pi/16;
                else
                    hi=pi/32;
                end%if
            end %if
        end %if

        Lsparsa=0; nzL=0; % per default lo Jacobiano si considera pieno

        % galleria di problemi-test

        if pb==1

            % varie possibilita' per Jacobiana scalare

            %IVPfun=1; % nome funzione: g1.m
            IVPfun=0;

            %L=-1e+5;
            %L=-1e+3

            %L=-4.804623167696255e+001, disp('max abs( eig(L)), pb==12, N=10')
            %L=-9.981363861300441e-001, disp(' min abs(eig(L)), pb==12, N=20')
            %L=-9.932212059289505e-001, disp('min abs(eig(L)), pb==12, N=10')
            %L=-1.777324315589538e+002, disp(' max abs( eig(L)), pb==12, N=20')

            %L=+9.836316430834666e-001i, disp('min(abs(eig(L)) per pb==11, N=10')
            %L=+3.183098861837908e+000i, disp('max abs( eig(L)), pb==11, N=10')

            %L=9.958927352435616e-001i, disp('min(abs(eig(L)) per pb==11, N=20')
            %=2*pi+00i, disp('max abs( eig(L)), pb==11, N=20')


            %L=-3.997664453664651e+000, disp('Scalar Jacobian matrix - large eig.')
            %L=2.335546335346055e-003, disp('Scalar Jacobian matrix - small eig.')


            L=input('scalar J=');
            %L=-1

            y0=1;
            t0=0; T=1;
        end


        if pb==2                          % oscillatore armonico
            %L=[0 1; -pi*pi 0];
            %L=[0 -1; 1 0]
            %gamm=input('gamma per osc. arm.=');
            omeg=1; omeg2=omeg*omeg;
            L=[0 1; -omeg2 0]; y0=[1 0];
            t0=0; T=2*pi;
            %T=input('Tempo finale T=');
        end

        if pb==3                        % L. Lopez, D. Trigiante
            L=[-21 19 -20;                % appl. num. math. 11 (1993)
                19 -21 20;                % 225-239
                40 -40 -40];
            y0=[1 0 -1];
            t0=0;
            T =2.5;
        end

        if pb==4
            IVPfun=1; disp('Usare GBDF, pb stiff'); % nome funzione: g4.m
            L=[   -2    1  ;
                998 -999 ];               % J. D. Lambert,
            y0=[2 3];                      % Num.  Methods for ODEs (1991)
            t0=0; T=10;
        end

        %if pb==5
        %  L=[ -150.5 -149.5;   % sbagliato - vedi articolo
        %       149.5 -150.5 ];               % Lopez, Trigiante
        %  y0=[0 2];                          % appl. num. math. 11 (1993)
        %  t0=0; T=1;
        %  disp('sbagliato - vedi articolo');
        %end
        if pb==5
            disp('yp=y-2*sin(t), y0=1')
            IVPfun=5;
            L=1;
            y0=1;
            t0=0; T=input('Final T=');
        end

        if pb==22 % ACCELEROMETRO PROPOSTO DAGLI INGEGNERI DI FIORELLA
            IVPfun=22; gg=9; dd=1; % parametri dell'eq
            L=[0 1; gg/dd 0]; y0=[1 0];
            t0=0; T=60; %secondi
            %T=input('Tempo finale T=');
        end

        if pb==6                                % L. Brugnano, D. Trigiante
            L=-[1e4 1e3 1e2; 0 1e1 1e0; 0 0 1e-5]; % appl. num. math. 13 (1993)
            y0=[1 0 -1]; % initial condition       % pp. 277-290
            t0=0;
            T=2.50;
        end

        if pb==7                    % J. D. Lambert
            L=[42.2  50.1  -42.1;     % Num. Methods for ODEs, pp. 228-229
                -66.1 -58.0   58.1;     % (modi debolm. crescenti)
                26.1  42.1  -34.0 ];
            y0=[1  0  2];
            t0=0; T=3;
        end

        if pb==9                     % problema Hamiltoniano full 10x10
            for i=1:10                  % L=[0 -I_5; I_5 0]*8 I_10 + Q
                for j=1:10                % P. Amodio, L. Brugnano
                    Q(i,j)=i+j;             % JCAM 78 (1997) 197-211
                end
            end
            L=sparse([zeros(5) -eye(5); eye(5) zeros(5)])*(8*speye(10)+Q);
            y0=[1 zeros(1, 9)];
            t0=0; T=20;
        end

        if pb==10
            titolo=['advenction equation u_t+c*u_x=0'];
            N=input('Introduci N:');
            x_max=pi; t0=0; c=1;         % advenction equation u_t+c*u_x=0
            T=2*pi; %N=50;               % discretizzata con forward Euler
            dx=x_max/(N+1); %dt=h=T/s    % per LMOL (longitudinal m.o.l.)
            y0=zeros(1,N);               % u(x,0)=g(x)=sin(x)
            L=spalloc(N,N,2*N-1);        % u(x_max,t)=0
            for i=1:N                    % Dominio=[0,x_max]x[0,T]
                y0(i)=i*dx; % u(x,0)=x
                %y0(i)=dx*sqrt(i*(N-i));
                L(i,i)=-c/dx;
                if (i<N), L(i,i+1)=c/dx; end
            end %i
            %y0(N/2)=1e+6; % picco nella condizione iniziale!!
            fprintf('\nAdvection eq., bc=0, mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; nzL=2*N; nzdiag=2;
        end


        if pb==11
            titolo='advenction equation (diff. centr.)';
            N=input('Introduci N:');     % advenction equation u_t+c*u_x=0
            x_max=3; t0=0; c=1;         % discretizzata con differenze centrate.
            T=6;                      % Condizioni al contorno periodiche
            dx=x_max/(N); %dt=h=T/s    % (longitudinal m.o.l.)
            y0=zeros(1,N);               % u(x,0)=g(x)=sin(x)
            L=[]; L=spalloc(N,N,2*N);    % u(x_max,t)=u(0,t)=?
            dx2=2*dx;                    % Dominio=[0,x_max]x[0,T]
            for i=1:N

                %y0(i)=sin((i-1)*dx);  % vecchie cond. iniziali
                y0(i)=dx*dx*i*(N-i); % u(x,0)=x(x_{max}-x)

                if (i>1), L(i,i-1)=c/dx2; end;
                if (i<N), L(i,i+1)=- c/dx2; end;
            end %i
            L(1,N)=-c/dx2; % ci vuole il '-'
            L(N,1)=c/dx2;  % se messo=0, allora  caso B.C. non periodiche

            fprintf('\nAdvection eq. periodic b.c. mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; %onceLU=0
        end


        if pb==12
            titolo='equazione del calore u_t-u_{xx}=0';
            N=input('Introduci N:');
            x_max=pi; t0=0;              % equazione del calore u_t-u_{xx}=0
            T=2*pi;                      % discretizzata con LMOL
            dx=x_max/(N+1); %dt=h=T/s    % (longitudinal m.o.l.)
            y0=zeros(1,N);               % u(x,0)=g(x)=sin(x)
            for i=1:N                    % u(0,t)=u(x_{max},t)=0, 0<=t<=T
                %y0(i)=sin(i*dx);           % Dominio=[0,x_max]x[0,T]
                %y0(i)=rand;
                y0(i)=i*dx; % u(x,0)=x
            end %i
            %y0(N/2)=1e+6; % picco nella condizione iniziale!!
            L=(1/(dx*dx))*sparse(toeplitz([-2 1 zeros(1,N-2)]));
            fprintf('\nEq. calore su mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; nzL=2*N;
        end

        if pb==13
            titolo='Eq. calore 2 dim.';
            N=input('Introduci N:');
            x_max=pi; t0=0;         % equazione del calore u_t=u_{xx}+u_{yy}
            T=2*pi; dx=x_max/(N+1); % bidim. in dominio rettangolare (LMOL)
            d2=1/(dx*dx); N2=N*N;   % \omega:=[0,x]x[0,y]x[0,T]
            y0=zeros(1,N2);         % u((x,y),0)=g(x,y)=sin(x)*sin(y)
            L=spalloc(N2,N2,5*N2-4*N);% u((x,y),t)=0, 0<=t<=T, (x,y)\in \omega
            for j1=0:N-1
                for j2=1:N
                    rr=j1*N+j2;
                    %d2=d2/abs((j1+1)*dx+j2*dx); % fattore smorzante
                    L(rr,rr)=-4*d2;
                    if (j2>1), L(rr,rr-1)=d2; end
                    if (j2<N), L(rr,rr+1)=d2; end
                    if (j1<N-1), L(rr,rr+N)=d2; end
                    if (j1>0), L(rr,rr-N)=d2; end
                    %y0(rr)=((-1)^(j1))*(j1^3*tan(j1*dx)*j2);
                    y0(rr)=(j1*dx*j2*dx); % condiz. iniziale
                    %y0(rr)=rand;
                end %j2
            end %j1
            %y0(N2/2)=1e+6;
            fprintf('\nEq. calore bidimensionale, mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; nzL=2*N;
        end


        if pb==14
            titolo='Equazione delle onde u_{tt}-u_{xx}=0';
            N=input('Introduci N:');     % Equazione delle onde u_{tt}-u_{xx}=0
            x_max=pi; t0=0; T=2*pi;      % discretizzata con LMOL
            dx=x_max/(N+1); %dt=T/s      % (longitudinal m.o.l.)
            Np=N+1; dx2=dx*dx; N2=N*N;
            L=spalloc(2*N,2*N,4*N-2);
            y0=zeros(1,2*N);             % u(x,0)=g(x)=x, u_t(x,0)=0, 0<=x<=x_{max}
            for i=1:N                    % u(0,t)=u(x_{max},t)=0, 0<=t<=T
                % rappr. nel dominio=[dx,x_{max}]x[0,T]
                %y0(i)=(i*(N-i)*dx)/N;
                %y0(i)=0; % u(x,0)=0
                y0(i+N)=i*dx; % u(x,0)=x (da' sol. stabile)
                %y0(i)=i*dx; if (i==1), disp('ATX! u(x,0)=x'), end

                L(i,N+i)=1;
                L(N+i,i)=-2/dx2;
                if (i>1), L(N+i,i-1)=1/dx2; end %if
                if (i<N), L(N+i,i+1)=1/dx2; end %if
            end %i
            %y0(N/2)=1; % picco nella c.i.
            fprintf('\nEq. ONDE su mesh %i x %i (N x s+1)',N,s)
            Lsparsa=1; nzL=N*N; nzdiag=3;
        end

        if pb==15
            titolo='advenction equation (diff. centr. 2)';
            N=input('Introduci N:');     % advenction equation u_t+c*u_x=0
            x_max=pi; t0=0; c=1;         % discretizzata con differenze centrate.
            T=2*pi;                      % Condizioni al contorno periodiche
            dx=x_max/(N); %dt=h=T/s    % (longitudinal m.o.l.)
            y0=zeros(1,N);               % u(x,0)=g(x)=sin(x)
            L=[]; L=spalloc(N,N,2*N);    % u(x_max,t)=u(0,t)=?
            dx2=2*dx;                    % Dominio=[0,x_max]x[0,T]
            for i=1:N
                %y0(i)=sin((i-1)*dx);
                y0(i)=(i*(N-i)*dx);
                if (i>1), L(i,i-1)=-c/dx2; end;
                if (i<N), L(i,i+1)= c/dx2; end;
            end %i
            L(1,N)=c/dx2;
            L(N,1)=-c/dx2;  % se messo=0, allora
            nzdiag=5;                               % caso B.C. non periodiche
        end

        if pb==16
            titolo='eq. calore con coeff. var.';
            N=input('Introduci N:');
            x_max=pi; t0=0;   T=2*pi;
            n=3; c1=1; c2=0; %coefficienti per v.i., jacob.
            dx=x_max/(N+1);
            y0=zeros(1,N);
            L=spalloc(N,N,3*N);
            for i=1:N
                %y0(i)=sin(c1*dx*i+c2); if (i==1), disp('g(x)=sin(x)'), end
                y0(i)=(dx*i); if (i==1), disp('g(x)=x'), end
                dx1=exp(-(dx*(i-.5))^n);
                dx2=exp(-(dx*(i+.5))^n);
                if (i>1), L(i,i-1)=dx1; end;
                if (i<N), L(i,i+1)=dx2; end;
                L(i,i)=-dx1-dx2;
            end %i
            L=(1/(dx*dx))*sparse(L);
            fprintf('\nEq. calore su mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; nzL=2*N; nzdiag=3;
        end

        if pb==18 % eq. diffusione 2D diff. centrate con coeff. diffusione var.
            titolo='Eq. diff. 2D coeff-var.';
            N=input('Introduci N:');
            x_max=3; t0=0;         % equazione del calore u_t=u_{xx}+u_{yy}
            T=6; dx=x_max/(N+1); % bidim. in dominio rettangolare
            d2=1/(dx*dx); N2=N*N;   % \omega:=[0,x]x[0,y]x[0,T]
            y0=zeros(1,N2);         % u((x,y),0)=g(x,y)
            L=spalloc(N2,N2,5*N2-4*N);% u((x,y),t)=0, 0<=t<=T, (x,y)\in \omega

            %d2=(50+1)^2/((x_max)^2), disp('ATX! SOLO PER GRAFICO')

            for j1=0:N-1
                for j2=1:N
                    rr=j1*N+j2; % natural ordering
                    j1p=j1+1; ru=j1p*dx;
                    aa1=aa((j2-.5)*dx,j1p*dx)*d2; aa2=aa((j2+.5)*dx,j1p*dx)*d2;
                    aa3=aa(j2*dx,(j1p+.5)*dx)*d2; aa4=aa(j2*dx,(j1p-.5)*dx)*d2;
                    L(rr,rr)=-aa1-aa2-aa3-aa4;
                    if (j2>1),   L(rr,rr-1)=aa1; end
                    if (j2<N),   L(rr,rr+1)=aa2; end
                    if (j1<N-1), L(rr,rr+N)=aa3; end
                    if (j1>0),   L(rr,rr-N)=aa4; end
                    %y0(rr)=((-1)^(j1))*(j1^3*tan(j1*dx)*j2);
                    y0(rr)=(j1*dx*j2*dx); % condiz. iniziale
                    %y0(rr)=rand;
                end %j2
            end %j1
            fprintf('\nEq diff. 2D coeff var, mesh %i x %i (N x s)',N,s)
            Lsparsa=1; nzL=2*N;
        end

        if pb==19
            a=input('Re(l):');
            b=input('Im(l):');
            L=[a b; -b a]; y0=[1 0];
            t0=0; T=2*pi;
            %T=input('Tempo finale T=');
        end

        if pb==20
            titolo='Laplaciano puro';
            N=input('Introduci N:');
            x_max=1; t0=0;              % equazione del calore u_t-u_{xx}=0
            T=1;                      % discretizzata con LMOL
            dx=x_max/(N+1); %dt=h=T/s    % (longitudinal m.o.l.)
            y0=zeros(1,N);               % u(x,0)=g(x)=sin(x)
            for i=1:N                    % u(0,t)=u(x_{max},t)=0, 0<=t<=T
                %y0(i)=sin(i*dx);           % Dominio=[0,x_max]x[0,T]
                %y0(i)=rand;
                y0(i)=i*dx; % u(x,0)=x
            end %i
            %y0(N/2)=1e+6; % picco nella condizione iniziale!!
            L=sparse(toeplitz([-2 1 zeros(1,N-2)]));
            fprintf('\nEq. calore su mesh %i x %i (N x s+1)',N,s+1)
            Lsparsa=1; nzL=2*N;
        end

        %% PROBLEMI FRAZIONARI AGGIUNTI DA FABIO:

        if pb==21   % Problema: Example 1 da Articolo Cinesi [Gu et Al.]
            addpath('../generatori_FDE/')
            titolo='Fractional Diffusion Equation';
            %N=input('Introduci N:');
            N = Nvalue;
            x0 = 0; y0 = 0;
            x_max = 2; yN = 0;
            l = @(x) 0.6*ones(size(x));
            r = @(x) 0.5*ones(size(x));
            b = @(x) zeros(size(x));
            c = @(x) zeros(size(x));
            fracexp = 1.2;
            T = 1;
            t0 = 0;
            [L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);
            L = -L;
            xx = x0:dx:x_max;
            g21 = @(t) (-xx(2:end-1).*(xx(2:end-1) - 1))' + F;
            IVPfun = 0;
            y0 = exp(-(xx(2:end-1) - 1.2).^2./(2*(0.08)^2));
            % Definizione del precondizionatore per lo Jacobiano (LowMemory)
        end

        if pb==22   % Problema: Example 1 da Articolo Cinesi con altro alpha [Gu et Al.]
            addpath('../generatori_FDE/')
            titolo='Fractional Diffusion Equation';
            %N=input('Introduci N:');
            N = Nvalue;
            x0 = 0; y0 = 0;
            x_max = 2; yN = 0;
            l = @(x) 0.6*ones(size(x));
            r = @(x) 0.5*ones(size(x));
            b = @(x) zeros(size(x));
            c = @(x) zeros(size(x));
            fracexp = 1.5;
            T = 1;
            t0 = 0;
            [L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);
            L = -L;
            xx = x0:dx:x_max;
            IVPfun = 0;
            g22 = @(t) (-xx(2:end-1).*(xx(2:end-1) - 1))' + F;
            y0 = exp(-(xx(2:end-1) - 1.2).^2./(2*(0.08)^2));
            % Definizione del precondizionatore per lo Jacobiano (LowMemory)
        end

        if pb==23   % Secondo problema di Test dei Cinesi [Gu et Al.]
            addpath('../generatori_FDE/')
            titolo='Fractional Diffusion Equation';
            %N=input('Introduci N:');
            N = Nvalue;
            %fracexp = input('Introduci alpha := ');
            fracexp = fracval;
            x0 = 0; y0 = 0;
            x_max = 2; yN = 0;
            l = @(x) eulgamma(3-fracexp).*x.^fracexp;
            r = @(x) eulgamma(3-fracexp).*(2-x).^fracexp;
            b = @(x) zeros(size(x));
            c = @(x) zeros(size(x));
            T = 1;
            t0 = 0;
            [L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);
            L = -L;
            xx = x0:dx:x_max;
            IVPfun = 1;
            g23 = @(t) (-32*exp(-t).*(...
                xx(2:end-1).^2 + 1/8.*(2-xx(2:end-1)).^2.*(8+xx(2:end-1).^2) ...
                -3/(3-fracexp).*(xx(2:end-1).^3 + (2-xx(2:end-1)).^3) + 3/((4-fracexp)*(3-fracexp)).*...
                (xx(2:end-1).^4 + (2-xx(2:end-1)).^4))+ F')';
            y0 = 4.*xx(2:end-1).^2.*(2-xx(2:end-1)).^2;
            % Definizione del precondizionatore per lo Jacobiano (LowMemory o BCCB-Strang)
        end

        if pb==24   % Problema: Example 1 da Articolo Cinesi con altro alpha [Gu et Al.]
            addpath('../generatori_FDE/')
            titolo='Fractional Advection Dispersion Equation';
            %N=input('Introduci N:');
            N = Nvalue;
            x0 = 0; y0 = 0;
            x_max = pi; yN = 0;
            l = @(x) 3*(1 + x.^2);
            r = @(x) 3*(1 + x.^2);
            b = @(x) 1+0.5*sin(4*x);
            c = @(x) 7-4/3*x;
            fracexp = fracval;
            T = 1;
            t0 = 0;
            [L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);
            L = -L;
            xx = x0:dx:x_max;
            IVPfun = 1;
            g24 = @(t) exp(-t).*sin(xx(2:end-1)') - F;
            y0 = -xx(2:end-1).*(xx(2:end-1) - pi);
            % Definizione del precondizionatore per lo Jacobiano (LowMemory)
        end

        if pb==25   % Problema: Example 1 da Articolo Cinesi con altro alpha [Gu et Al.]
            addpath('../generatori_FDE/')
            addpath('../MFC2/')
            titolo='Fractional Advection Dispersion Equation';
            N = Nvalue;
            T = 1; t0 = 0;
            a = 0;      % Physical domain [a,b]x[c,d]
            b = 1;
            c = 0;
            d = 1;
            alpha = 1.1;    % Fractional order(s) of derivative
            beta = 1.8;
            Kx1 = 2;   % Fractional diffusion coefficients
            Kx2 = 1.5;
            b1 = @(x,y) beta+0.5*sin(4*pi*x).*cos(5*pi*y);
            b2 = @(x,y) alpha+0.7*sin(7*pi*y).*cos(4*pi*x);
            cc = @(x,y) 1+0.5*cos(x.*y);
            % Discretization parameter:
            Nx1 =   N;     % Number of grid point in the x1-direction
            Nx2 =   N;     % Number of grid point in the x2-direction

            % Building of the discretization:
            x1 = linspace(a,b,Nx1); hx1 = x1(2)-x1(1);
            x2 = linspace(c,d,Nx2); hx2 = x2(2)-x2(1);
            [X1,X2] = meshgrid(x1(2:end-1),x2(2:end-1));
            ex1 = ones(Nx1-2,1);
            ex2 = ones(Nx2-2,1);
            Tx1 = -Kx1*ranort(alpha,Nx1,hx1);
            Tx2 = -Kx2*ranort(beta,Nx2,hx2);
            Ix1 = speye(Nx1-2);
            Ix2 = speye(Nx2-2);
            Bx1 = spdiags([-ex1./(2*hx1),ex1./(2*hx1)],[-1,1],Nx1-2,Nx1-2);
            Bx2 = spdiags([-ex2./(2*hx1),ex2./(2*hx2)],[-1,1],Nx2-2,Nx2-2);
            Dx1 = spdiags(reshape(b1(X1,X2),(Nx1-2)*(Nx2-2),1),0,(Nx1-2)*(Nx2-2),(Nx1-2)*(Nx2-2));
            Dx2 = spdiags(reshape(b2(X1,X2),(Nx1-2)*(Nx2-2),1),0,(Nx1-2)*(Nx2-2),(Nx1-2)*(Nx2-2));
            C = cc(X1,X2);

            L = -(kron(Tx1(2:end-1,2:end-1),Ix2)+kron(Ix1,Tx2(2:end-1,2:end-1)) ...
                + Dx1*kron(Bx1,Ix2) + Dx2*kron(Ix1,Bx2) ...
                + spdiags(reshape(C,(Nx1-2)*(Nx2-2),1),0,(Nx1-2)*(Nx2-2),(Nx1-2)*(Nx2-2)));

            Zd_state = @(x,y) sin(pi*5*x).*sin(pi*5*y);
            Y0_state = (X1-1).*(X2-1).*X1.*X2;
            z = reshape(Zd_state(X1,X2),(Nx1-2)*(Nx2-2),1);
            IVPfun = 1;
            g25 = @(t) z.*exp(-t);
            y0 = reshape(Y0_state,1,(Nx1-2)*(Nx2-2));
        end

        if pb == 26
            % Space-variant fractional problem
            N = Nvalue;
            alpha = 1.5; % Order of the fractional derivative
            u0f = @(x) 5*x.*(2-x);
            hN = 2/(N-1); x = 0:hN:2;
            u0 = u0f(x).';
            dplus = @(x,t) gamma(3-alpha).*x.^alpha;
            dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
            G = glmatrix(N,alpha);
            Gr = G; Grt = G.';
            Dplus = diag(dplus(x,0));
            Dminus = diag(dminus(x,0));
            I = eye(N,N);
            e = ones(N,1);
            L = (Dplus*Gr + Dminus*Grt)/hN^alpha;
            T = 1;
            t0 = 0;
            IVPfun = 1;
            g26 = @(t) -x.*(x-2).*exp(-t);
            y0 = u0f(x);
        end

        %% fine galleria di problemi-test %%


        m=max(size(L));
        color=['y' 'm' 'r' 'c' 'g' 'b' 'w' 'yo' 'mo' 'co' 'ro' 'go' 'bo'];

        % passo variabile 'artificiale'
        %p_var=input('Passo variabile (1/0)');
        p_var=0;
        if (p_var==1)
            p=input('Numero blocchi (p):');
            if (p==16)
                r=1.9; % per p=16, pb=3 va bene
            elseif (p==32)
                r=1.34; % per p=32
            elseif (p==64)
                r=1.14; % per p=64
            elseif (p==8)
                r=4; % per p=8
            end
            %  r=1.062, p=128
            %   r=4,p=8 % p=8
            %   r=35, p=4 % p=4
            %   r=500, p=2,
            if pb==1
                h_tilde(1)=1e-2;
            else
                h_tilde(1)=1e-4;
            end
            T=t0;
            for i=2:p
                h_tilde(i)=h_tilde(i-1)*r;
                T=T+h_tilde(i);
            end
            disp('Mesh e punto finale:')
            h_tilde;
            T;
        else % passo costante
            p=1;
            %p=input('Numero blocchi (p)=');
            passo=(T-t0)/p;
            for i=1:p
                h_tilde(i)=passo;
            end
        end

        hi=passo; hi=h_tilde(1)/s; % PASSO COSTANTE
        %(e' il passo h dell'articolo su SISC 22-3, 2000)
        if (hi<=0)
            hi=input('h= [se <0, h=(T-t0)/s]');
        end %if
        T=t0+hi*s; % da qui si puo' immettere h dall'esterno
        if hi<0
            hi=(T-t0)/s;
        end

        %fprintf('\nPb=%i tol=%e \n[%e,%e] precond=%i',pb,tol,t0,T,precond)
        %fprintf('\n step size=%e, residuo:%i, metodo:%i\n',hi,residuo_vero,metodo);
        % MAX NUMERO DI ITERAZIONI AMMESSE
        max_it=1000;

        fi=2; plot_error=0;
        msp=m*(s*p+1); % dimensione sistema a blocchi
        sp=s*p;

        %
        % creazione della matrice sparsa a blocchi Ax=b
        % del BVM opportuno
        %

        % ALPHA,BETA matrici dei coeff. del metodo, sono s x (s+1)
        [ALPHA,BETA]=mab(metodo, passi, s);

        A1=ALPHA(:,2:s+1); % costruisce le matrici quadrate che corrispondono
        B1=BETA(:,2:s+1);  % ad A e B come sull'articolo (queste sono s x s)

        if metodo==0
            k1 = fix( (passi+2)/2 );
        else
            k1 = ceil( passi/2 );
        end

        if metodo==3 % coefficienti metodo LMF determinati a parte per il precondizionatore)
            [ro,si] = rosi_tom( passi ); %disp('TOM')
        else
            [ro,si] = coef( metodo, passi, k1 );
        end % if

        if (metodo==-1 && passi==2) % Midpoint (o Simpson) +Implicit Euler di Axelsson
            simpson=input('Simpson (1) o Midpoint(0):'); % se posto a 1 da' Simpson invece di Midpoint
            ALPHA=2*ALPHA;
            ALPHA(s,s-1)=0; ALPHA(s,s)=-1; ALPHA(s,s+1)=1;
            if simpson
                for j=1:s-1, BETA(j,j)=1/3; BETA(j,j+1)=4/3; BETA(j,j+2)=1/3; end
                %BETA(s,s)=1/2; BETA(s,s+1)=1/2;
                BETA(s,s)=0; BETA(s,s+1)=1; % for Implicit Euler
                si(1)=1/3; si(2)=4/3; si(3)=1/3;
                %disp('ATX: Simpson+Implicit Euler.')
            else
                for j=1:s-1, BETA(j,j+1)=2; end
                BETA(s,s+1)=1; si(2)=2;
                %disp('ATX: MP+IE attivato.')
            end %if_simpson
            ro(1)=-1; ro(3)=1;
        end %if

        if (metodo==1 && passi==2) % Adams+Trapezi di Axelsson
            %[ALPHA BETA]=mab(1, 2, s);

            A1=ALPHA(:,2:s+1); % costruisce le matrici quadrate che corrispondono
            B1=BETA(:,2:s+1);  % ad A e B come sull'articolo (queste sono s x s)
            % si cambiano solo le condizioni finali

            BETA(s,s)=0; BETA(s,s+1)=1; %disp('ATX: Adams+Impl. Euler'), % for Implicit Euler
            BETA(s,s-1)=0;
            %BETA(s,s)=1/2; BETA(s,s+1)=1/2; disp('ATX: Adams+Trapezi'); % for T5rapezoidal rule
            %disp('ATX: Adams+Impl. Euler')
        end %if

        % MEMORY ALLOCATION
        %         if (pb<10) % vecchio tipo
        %           A=spalloc(msp,msp,m+p*m^2*s^2+p*m^2*s);
        %         else
        % %         fprintf('\nATTENZIONE: MAX % diagonali non zero nello JACOBIANO!\n',nzdiag);
        %          A=spalloc(msp,msp,m+nzdiag*(passi+1)*m*s+nzdiag*m);
        %         end %if
        %                 A=spalloc(msp,msp,(1+p*(passi+1)*s)*nnz(L));
        %
        %                 for i=1:m
        %                     A(i,i)=1; % Im nel primo blocco diag.
        %                 end
        %
        %                 if p==1
        %                     A=kron(ALPHA,speye(m))-hi*kron(BETA,L);
        %                 else
        %                     % creazione blocchi matrice di iterazione M
        %                     % (primo macroblocco)
        %                     ViMi=kron(ALPHA,speye(m))-hi*kron(BETA,L);
        %                     for i=m+1:m+m*s
        %                         for j=1:m+m*s
        %                             A(i,j)=ViMi(i-m,j); % primo blocco utile
        %                         end % j
        %                     end %i
        %                     % creazione successivi macroblocchi per block-BVM
        %                     for bl=2:p % dal blocco 2 a p
        %                         hi=h_tilde(bl)/s;
        %                         ViMi=kron(ALPHA,speye(m))-hi*kron(BETA,L);
        %                         r_o=m+1+(bl-1)*m*s;
        %                         c_o=(bl-1)*m*s+1;
        %                         for i=r_o:r_o+m*s-1
        %                             for j=c_o:c_o+m*s+m-1
        %                                 A(i,j)=ViMi(i-r_o+1,j-c_o+1); % blocchi da 2 a p
        %                             end % j
        %                         end %i
        %                     end %bl
        %
        %                 end % if p==1
        %        A = kron(ALPHA,speye(m))-hi*kron(BETA,L); %% BUILD ONLY TO TEST
        %        DIRECT SOLUTION
        Am = kron(ALPHA(:,1),speye(m))-hi*kron(BETA(:,1),L);
        %fprintf('\nK1*(A)=%e,K1*(B)=%e\n',condest((ALPHA(:,2:s+1))),condest((BETA(:,2:s+1))));
        %fprintf('\nK1*(M)=%e\n',condest((A(:,m+1:msp))));
        %disp('PREMI UN TASTO PER RISOLVERE IL SITEMA LIN.')
        %pause

        %%  SETUP PER VETTORI X e b DEL SISTEMA LINEARE
        %clear ALPHA BETA

        b=zeros(msp,1);
        % b=spalloc(msp,1,m); % se f(t)=0 nell'IVP, usare questa!!!

        if (IVPfun>0)
            pbstr=int2str(pb);
            for j=1:s
                b(1+j*m:(j+1)*m)=eval(['g' pbstr '(t0+j*hi)']);
            end %for j
            b(m+1:m*(s+1))=hi*kron(BETA,speye(m))*b;
        end %if
        b(1:m)=y0;

        %fprintf(1,'\nMatrice M (%i x %i)\n\n',m*s*p,m*s*p);

        %
        % PRECONDIZIONAMENTO DELLA MATRICE DI ITERAZIONE
        %


        if (precond==1) %
            % precondizionamento e stima di x(0) con Implicit Euler
            % (Brugnano-Mazzia-Trigiante 1996)

            x(1:m)=y0;
            P_tilde=spalloc(msp,msp,m+m*m*s*p);
            %P_tilde(1:m,1:m)=speye(m);
            P_tilde(1:m,1:m)=inv( eye(m)-((h_tilde(1)/s)*L) );
            for i=1:sp
                index=m*(i-1)+1;
                index2=m*i;
                j=ceil(i/s);
                I_hL_1=inv( eye(m)-((h_tilde(j)/s)*L) ); %(I-hi*L)^-1
                index3=m+1+m*(i-1);
                index4=index3+m-1;
                P_tilde(index3:index4,index3:index4)=I_hL_1;
                x(index2+1:index2+m)=I_hL_1*x(index:index2);
            end %i
            M=spalloc(msp,msp,m+m*m*s*p);
            P_1=spalloc(msp,msp,m+m*s*p+m*m*s*p);
            E=spalloc(msp,msp,m*s*p);
            P_1=speye(msp);
            for i=1:m*s*p
                E(m+i,i)=1;
            end %i
            P_1=P_tilde*E;           % Horner's algorithm
            M=P_1+speye(msp);
            for i=(p/2-2):-1:0
                M=P_1*M+speye(msp);
            end %i
            P_1=M*P_tilde;           % P_1 matrice di precond. invertita
            clear P_tilde M E I_hL_1
            A=P_1*A; % precondizionamento di A a sx
            b=P_1*b;
            %disp('Precondizionamento con BACKWARD EULERO avvenuto (1)')
            clear P_1

        elseif (precond==3 || precond==5)

            % PRECONDIZIONATORI STRANG-CIRCOLANTI A BLOCCHI
            % N.B.: s(A) E' SINGOLARE PER OGNI k
            % Bertaccini-Strang, 1997
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$ for T.Chan-circ
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            for j=0:s-1
                jp=j+1; i1=jp+k1; i2=jp-s+k1; aj_ro=0; aj_si=0;
                if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp), aj_ro=ro(i2); aj_si=si(i2); end
                r1(jp)=aj_ro; r2(jp)=aj_si;
            end %for_j
            if precond==3
                %disp('**STRANG circulant (3) **')
            else
                %disp('**MODIFIED STRANG circulant (5) **')
            end %if

        elseif (precond==6)  % ******
            % PRECONDIZIONATORI skew-CIRCOLANTI
            % (c) Bertaccini-Ng, 2000
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            omega=-1; % per caso skew-circolante
            %omega=input('Omega (circulant)=');
            OMEGA=omega.^(-(0:s-1)/s);
            for j=0:s-1
                jp=j+1; i1=jp+k1; i2=jp-s+k1; aj_ro=0; aj_si=0;
                if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp), aj_ro=-ro(i2); aj_si=-si(i2); end
                r1(jp)=aj_ro; r2(jp)=aj_si;
            end %for_j
            %disp('** {skew}-circulant (6) **')

        elseif (precond==7)  % ****** ###
            % PRECONDIZIONATORI OMEGA-CIRCOLANTI A BLOCCHI
            % costruzione delle matrici componenti w(A), w(B), prima riga
            % (c) Bertaccini-Ng 2001
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            %omega=-1; % per caso skew-circolante
            omega=input('Omega-circulant [1=circ;-1=skew]=');
            OMEGA=omega.^(-(0:s-1)/s);
            for j=0:s-1
                jp=j+1; i1=jp+k1; i2=jp-s+k1; aj_ro=0; aj_si=0;
                if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp), aj_ro=-ro(i2); aj_si=-si(i2); end
                r1(jp)=aj_ro; r2(jp)=aj_si;
            end %for_j
            %disp('** {omega}-circulant (7) **')


        elseif (precond==2)
            %
            % PRECONDIZIONATORI CIRCOLANTI A BLOCCHI
            % (Bertaccini-T. Chan, 1997)
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            for j=0:s-1
                aj_ro=0; aj_si=0; aj_sro=0; aj_ssi=0;
                jp=j+1; i1=j+k1p; i2=j-s+k1p;
                if (i1<=kp) aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp) aj_sro=ro(i2); aj_ssi=si(i2); end
                r1(jp)=( (s-j)*aj_ro+j*aj_sro )/s; % riga per la matrice \breve A_s
                r2(jp)=( (s-j)*aj_si+j*aj_ssi )/s; % riga per la matrice \breve B_s
            end %for_j
            %disp('**T.CHAN circulant (type=2) **')

        elseif (precond==4)
            % NUOVA CLASSE DI PRECONDIZIONATORI CIRCOLANTI A BLOCCHI
            % (c) D. Bertaccini, 1997
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$ for DAN-CIRC
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            for j=0:s-1
                aj_ro=0; aj_si=0; aj_sro=0; aj_ssi=0;
                i1=j+k1p; i2=j-s+k1p; jp=j+1;
                if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp), aj_sro=ro(i2); aj_ssi=si(i2); end
                r1(jp)=( (s+j)*aj_ro+j*aj_sro )/s; % riga per la matrice \breve A_s
                r2(jp)=( (s+j)*aj_si+j*aj_ssi )/s; % riga per la matrice \breve B_s
            end %for_j
            disp('** P-circulants (type=4) **')

        elseif (precond==8) % Precondizionatore 'Felix'
            %
            % Coefficients of $\breve A_s, B_s$ for $\phi, \psi$ for T.Chan-circ
            gamma=1; % COEFFICIENTE PER QUESTA APPROX. CIRCULANTE
            sp=s+1; mp=m+1; ms=m*s; m2=m*2; msp=m*s*p;
            kp=passi+1; k1p=k1+1;
            for j=0:s-1
                jp=j+1; i1=jp+k1; i2=jp-s+k1; aj_ro=0; aj_si=0;
                if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                if (i2>0 && i2<=kp), aj_ro=ro(i2); aj_si=si(i2); end
                r1(jp)=aj_ro; r2(jp)=aj_si;
            end %for_j
            r1(1)=r1(1)+gamma/s; r2(1)=r2(1)+gamma/s;
            disp('**Felix-STRANG circulant (8) **')
        end %if

        % set-up per precondizionatori circolanti Bertaccini
        if (precond==2 || precond==3 || precond==4 || precond==5 || precond==6 || precond==7 || precond==8)

            % calcolo gli autovalori di $\breve A_s$, $\breve B_s$
            if (precond==6 || precond==7)
                phi=(s)*ifft((OMEGA).*r1); % ### autovalori di $\breve A_s$
                psi=(s)*ifft((OMEGA).*r2); % autovalori di $\breve B_s$
            else
                phi=(s)*ifft(r1); % autovalori di $\breve A_s$
                psi=(s)*ifft(r2); % autovalori di $\breve B_s$
            end %if

            % approssimazione circolante MODIFIED STRANG per A
            if precond==5
                phi(1)=1/(s); %psi(1)=1-(1/(s));
                %phi(1)=1/(s); % ***

                %phi(1)=real(phi(s)); psi(1)=real(psi(s));

                r1=real(fft(phi/s)); % ricovero dei coeff. originali
                r2=real(fft(psi/s)); % ATTENZIONE! vero solo per Jacobiana
            end %if                % a valori REALI!!

            normL=norm(full(L));
            %coe=1/s; % coefficiente di scalatura per FFTs  in dbprec
            %coe=1/sqrt(s);
            coe=1;
            sqL=1e+10;

            % costruisce i vettori di permutazione e di antipermutazione
            index1=zeros(1,m*s*p); % permutazione
            index2=zeros(1,m*s*p); % antipermutazione
            i3=0;
            for i1=0:p-1              % permutazione
                for j=0:m-1
                    for i2=0:s-1
                        i3=i3+1;
                        index1(i3)=1+m*i2+j+i1*ms;
                    end
                end
            end
            i3=0;
            for i1=0:p-1              % antipermutazione
                for j=0:s-1
                    for i2=0:m-1
                        i3=i3+1;
                        index2(i3)=1+s*i2+j+i1*ms;
                    end
                end
            end

            % global normL phi psi coe index1 index2 h_tilde hi residuo_vero precond
            % global L sqL m p s Lsparsa onceLU nzL apprJ pol_flag
            %global L1 U1 Lg Ug Jc apprJ
            % if (precond==6)
            %    global OMEGA
            % end %if
        end % fine set-up per precondizionatori circolanti Bertaccini

        %% PRECONDIZIONATORE STRANG PER Jm
        if apprJ == -1 && fracprec == 1
            %CL = circ_approx(L,1);
            CL = [L(1:ceil(m/2),1);L(1,m-ceil(m/2)+1:-1:2)'];
            Jc = fft(full(CL(:,1)));
            if any(abs(Jc) < 1e-5)      % Modifica dei Cinesi
                phi(1)=real(phi(s));
            end
        end
        %% PRECONDIZIONATORE STRANG PER Jm - Versione 2D
        if apprJ == 0 && fracprec == 6
            nn=N-2;
            n2=(nn-1)/2;
            Lb = @(iii,jjj) L(nn*(iii-1)+1:nn*iii,nn*(jjj-1)+1:nn*jjj); % Comando per accedere al blocco (i,j)
            for iii=1:n2+1
                Atemp = Lb(iii,1);
                Stra=[Atemp(1:n2+1,1)' Atemp(1,nn-n2:-1:2)]; % Strang approximation - local
                col{iii} = gallery('circul',[Stra(1),Stra(end:-1:2)]);
            end
            for jjj=nn-n2:-1:2
                Atemp = Lb(1,jjj);
                Stra=[Atemp(1:n2+1,1)' Atemp(1,nn-n2:-1:2)]; % Strang approximation - local
                row{jjj} = gallery('circul',[Stra(1),Stra(end:-1:2)]);
            end

            Stra = {col{:},row{nn-n2:-1:2}};
            Stra = {Stra{1},Stra{end:-1:2}};
            NN = length(Stra);
            FPR = cell2mat(Stra(gallery('circul',(1:NN))));
        end
        %% PRECONDIZIONATORE SKEW-CIRCULANT PER Jm
        if apprJ == -1 && fracprec == 2
            apprJ = -1;
            CL = circ_approx(L,4);
            omega = -1;
            OMEGAm=omega.^(-(0:m-1)'/m);
            Jc = (m)*ifft(OMEGAm.*CL(:,1));
        end
        %% PRECONDIZIONATORE LOW-MEMORY PER Jm
        if apprJ == 0 && fracprec == 3
            if pb == 25
                diagN = ceil((N-2)^2/diagNdiv);
                FPR = spdiags( spdiags(L, -diagN:diagN), -diagN:diagN, (N-2)^2, (N-2)^2);
            else
                diagN = ceil(N/diagNdiv);
                FPR = spdiags( spdiags(L, -diagN:diagN), -diagN:diagN, N-1, N-1);
            end
        end
        if apprJ == 0 && (fracprec == 4 || fracprec == 5 || fracprec == 7)
            if pb == 25
                diagN = diagNdiv;
                FPR = spdiags( spdiags(L, -diagN:diagN), -diagN:diagN, (N-2)^2, (N-2)^2);
            else
                diagN = diagNdiv;
                FPR = spdiags( spdiags(L, -diagN:diagN), -diagN:diagN, N-1, N-1);
            end
        end
        % PER PLOTTARE SPETTRO PRECONDIZ./NON PRECONDIZ.
        %
        if (plot_spectra) % normalmente inattivo, solo per grafico spettro ++
            if precond==0 % per evitare assurdita'
                non_prec_plot=1;
            end %

            if (precond>0)
                if precond==6
                    r1c=r1; r2c=r2;
                    for j=0:s-1 % devo ricostruire il vettore r1 di strang per 'fare' c1
                        jp=j+1; i1=jp+k1; i2=jp-s+k1; aj_ro=0; aj_si=0;
                        if (i1<=kp), aj_ro=ro(i1); aj_si=si(i1); end
                        if (i2>0 && i2<=kp), aj_ro=ro(i2); aj_si=si(i2); end
                        r1(jp)=aj_ro; r2(jp)=aj_si;
                    end %for_j
                    iskew=ceil(s/2)+1;
                    c1=r1(2:s); c1(s)=r1(1); c1=c1(s:-1:1); % c prima colonna di OC1
                    c1(iskew:s)=-c1(iskew:s);
                    c2=r2(2:s); c2(s)=r2(1); c2=c2(s:-1:1); % c idem per OC2
                    c2(iskew:s)=-c2(iskew:s);
                    r1=r1c; r2=r2c;
                else
                    c1=r1(2:s); c1(s)=r1(1); c1=c1(s:-1:1); % c prima colonna di OC1
                    c2=r2(2:s); c2(s)=r2(1); c2=c2(s:-1:1); % c idem per OC2
                end %if precond==6
            end %if precond>0
            if (non_prec_plot && eig_computation)
                disp('Calcolo eigensystem (con precond) e condizionam.')
                [XX_nonp,DD]=eig(full(A(:,m+1:m*p*(s+1))));
                eigA=diag(DD);
            end %if
            if precond>0
                OC1=sparse(toeplitz(c1,r1)); % matrice per RO
                OC2=sparse(toeplitz(c2,r2)); % matrice per SIGMA

                if (pb==1 && 0) % normalmente disattivato.
                    OC=OC1/(-hi*L)+OC2; % matrice P nei pb scalari, forma alternativa
                    % per il solo calcolo della XX_p se hi*|L| e' grande
                else
                    OC=kron(OC1,speye(m))-kron(OC2,hi*L);% matrice P
                end %if
                %clear OC1 OC2

                % "precondizionamento" del calcolo degli autovalori
                % spesso migliora moltissimo il cond. dell'eigensystem
                % associato al cluster per skew e Strang
                % (quelli che... tight cluster)

                if (pb==1)
                    EE=(ALPHA(:,2:s+1)-OC1)-hi*L*(BETA(:,2:s+1)-OC2); % matrice di rango basso,
                    norOC1=norm(full(inv(OC2)*OC1))/(hi*abs(L));                 % problema scalare
                else
                    EE=A(:,m+1:m*p*(s+1))-OC; %
                end %if
                %PM=eye(m*s)+inv(OC)*EE;
                if (eig_computation)
                    disp('Calcolo eigensystem e suo numero cond...')
                    if ( (pb==1) && norOC1<1 )
                        PM=(OC1/(hi*L)+OC2)*inv(OC2)*EE;
                        disp('Approx al primo ordine di inv(P) applicata.')
                    else
                        PM=inv(OC)*EE; % da' un eigensystem meglio condizionato
                    end %if           % che con l'approccio "standard"
                    %DD=spalloc(msp,msp,msp);
                    [XX_p,DD]=eig(full(PM));
                    %[XX_p,DD]=eigs((PM));
                    eigPA=1+diag(DD);
                    clf; hold on; plot(eigPA,'r+')
                    %eigPAstd=eig(full(inv(OC)*A(:,m+1:m*p*(s+1))));
                    %plot(eigPAstd,'bo')
                    xlabel('Eigenv. of P^-1*M')
                    if (icn_yes)
                        YY_p=inv(XX_p);
                        for j=1:m*s*p % calcola i.c.n. per A
                            YY_p(j,:)=YY_p(j,:)/norm(YY_p(j,:)); y=YY_p(j,:);
                            %y=YY_p(j,:).'/norm(YY_p(j,:));
                            icn(j)=1/abs(y*XX_p(:,j));
                            fprintf('\n 1/icn(%i)=%e, l=(%16.15e,%e)',j,icn(j),real(eigPA(j)),imag(eigPA(j)));
                        end % for
                        fprintf('\n\nsum(1/icn(j))=%e\n',sum(icn));
                    end %if (icn_yes)
                    if non_prec_plot
                        fprintf('\nK2(V(M))=%e, \n',cond(XX_nonp));
                        subplot(1,2,1)
                        xlabel('Eigenv. of M')
                        plot(eigA,'g+')
                        fprintf('\n eigA: [%e,%e]x\nx[%e,%e]\n',min(real(eigA)),max(real(eigA)),min(imag(eigA)),max(imag(eigA)));
                        subplot(1,2,2)
                    end %if non_prec_plot
                    %hold on
                    hold off
                    fprintf('K2(V(P^-1*M))=%e\n',cond(XX_p));
                    fprintf('\n eigPA: [%e,%e]x\nx[%e,%e]\n',min(real(eigPA)),max(real(eigPA)),min(imag(eigPA)),max(imag(eigPA)));
                else
                    %fprintf('\nK1(OC)=%e, \n',condest(OC));
                    %filesuono='c:\windows\media\ding.wav';
                    %[Y,FS,BITS]=wavread(filesuono);
                    %sound(Y,FS,BITS);
                end %if eig_computation
            elseif non_prec_plot % non preconditioned case
                plot(eigA,'g+')
                xlabel('Eigenv. of M')
                fprintf('\nK2(V(M))=%e, \n',cond(XX_nonp));

            end %if precond>0

            %if precond==0
            %  fprintf('\nK1*(M)=%e\n',condest((A(:,m+1:msp))));
            %else
            %  fprintf('\nK1*(P)=%e\n',condest(OC));
            %end %if

            %eigPA=eig(full(inv(OC)*A(:,m+1:m*p*(s+1))));
            %eigA=eig(full(A(:,m+1:m*p*(s+1))));

            %svdPA=svd(full(inv(OC)*A(:,m+1:m*p*(s+1))));
            %svdA=svd(full(A(:,m+1:m*p*(s+1))));

            %filesuono='c:\windows\media\ding.wav';
            %[Y,FS,BITS]=wavread(filesuono);
            %sound(Y,FS,BITS);

        end %if (plot_spectra)

        %
        % SOLUZIONE SISTEMA LINEARE Mx=b
        %
        msp=m*(s*p+1); % dimensione sistema a blocchi
        %b(m+1:msp)=b(m+1:msp)-A(:,1:m)*y0';
        b(m+1:msp)=b(m+1:msp)-Am*y0';
        x=zeros(msp,1); %disp('x0=(0...0)');
        %x=ones(msp,1); disp('x0=(1...1)');

        %disp('SISTEMA PRONTO AL LANCIO...')
        %flops(0);
        addpath('../testfabio/');
        if (mspl==0) %BiCGSTAB(2)
            tic; % starts a stopwatch timer
            [x(m+1:msp), error, iter, iflag, res_iter0, m_v] = bicg2psx(@(x) Mprod(A1,B1,L,hi,x), x(m+1:msp), b(m+1:msp), max_it, tol,precond);
            tempo=toc; operaz=flops;
            if iflag ~= 0
                m_v = -m_v;
            end
        elseif (mspl==1) % GMRES(restart_it)
            if (precond<=0)
                tic; % starts a stopwatch timer
                % restart_it e' il numero di it. per ogni restart
                %             while (error>tolb && restarting<49)
                %                 [x(m+1:msp), error, iter, iflag, res_iter0] = gmrespsx( A(:,m+1:msp), x(m+1:msp), b(m+1:msp), restart_it, tol,precond);
                %                 restarting=restarting+1;
                %                 res_iter=[res_iter res_iter0];
                %             end %while
                %             iter=iter+restarting*restart_it;
                [x(m+1:msp),iflag,error,iter,res_iter0] = gmres(@(x) Mprod(A1,B1,L,hi,x), b(m+1:msp),restart_it,tol,max_it,[],[],x(m+1:msp));
                %fprintf('\n it. Restarted GMRES(%i)=%i',restart_it,iter)
                tempo=toc; operaz=flops; m_v=restart_it*iter(1)+iter(2);
                if iflag ~= 0
                    m_v = -m_v;
                end
            else
                %fprintf(1,'\n L.P. GMRES(%i) \n',restart_it)
                %if
                %restarting=-1; error=1;
                tolb=tol*norm(b(m+1:msp)); res_iter=[];
                tic; % starts a stopwatch timer
                % restart_it e' il numero di it. per ogni restart
                %             while (error>tolb && restarting<49)
                %                 [x(m+1:msp), error, iter, iflag, res_iter0] = gmrespsx( A(:,m+1:msp), x(m+1:msp), b(m+1:msp), restart_it, tol,precond);
                %                 restarting=restarting+1;
                %                 res_iter=[res_iter res_iter0];
                %             end %while
                %             iter=iter+restarting*restart_it;
                [x(m+1:msp),iflag,error,iter,res_iter0] = gmres(@(x) Mprod(A1,B1,L,hi,x), b(m+1:msp),restart_it,tol,max_it,@(v) dbprec(v),[],x(m+1:msp));
                %fprintf('\n it. Restarted GMRES(%i)=%i',restart_it,iter)
                tempo=toc; operaz=flops; m_v=restart_it*iter(1)+iter(2);
                if iflag ~= 0
                    m_v = -m_v;
                end
            end
        elseif (mspl==2) % (left preconditioned) BiCGStab
            if (precond<=0)
                %fprintf('\nunprecond. BiCGStab\n')
            else
                %fprintf('\nLeft precond. BiCGStab\n')
            end %if
            % starts a stopwatch timer
            %[x(m+1:msp), error, iter, iflag, res_iter, m_v]=bicgspsx(A(:,m+1:msp), x(m+1:msp), b(m+1:msp), max_it, tol,precond);
            if (precond > 0)
                %[x(m+1:msp),iflag,error,iter,res_iter] =  bicgstab(A(:,m+1:msp),b(m+1:msp),tol,max_it,@(v) dbprec(v),[],x(m+1:msp));
                tic;
                [x(m+1:msp),iflag,error,iter,res_iter] =  bicgstab(@(x) Mprod(A1,B1,L,hi,x),b(m+1:msp),tol,max_it,@(v) dbprec(v),[],x(m+1:msp));
                tempo=toc;
                if iflag ~= 0
                    iter = -iter;
                end
            else
                %[x(m+1:msp),iflag,error,iter,res_iter] =  bicgstab(A(:,m+1:msp),b(m+1:msp),tol,max_it,[],[],x(m+1:msp));
                tic;
                [x(m+1:msp),iflag,error,iter,res_iter] =  bicgstab(@(x) Mprod(A1,B1,L,hi,x),b(m+1:msp),tol,max_it,[],[],x(m+1:msp));
                tempo=toc;
                if iflag ~= 0
                    iter = -iter;
                end
            end
            operaz=flops; m_v = 2*iter;
            %fprintf('Iterazioni: %1.1f\n, Flag = %d\n',iter,iflag);
        elseif (mspl==3) % FGMRES+dbprec(v) with variable solver
            %fprintf(1,'\nBiCG...\n')
            tic; % starts a stopwatch timer
            %[x(m+1:msp), error, iter, iflag, res_iter,m_v]=bi_cg(A(:,m+1:msp), x(m+1:msp), b(m+1:msp), max_it, tol,precond);
            [x(m+1:msp),iflag,res,iter] = fgmres(@(x) Mprod(A1,B1,L,hi,x),b(m+1:msp),1000,max_it,tol);
            tempo=toc; operaz=flops;
            if iflag ~= 0
                iter = -iter;
            end
            m_v = iter;
        elseif (mspl==4) % BiCGStab(2)
            if (precond<=0)
                fprintf('\nunprecond. BiCGStab(2)\n')
            else
                fprintf('\nLeft precond. BiCGStab(2)\n')
            end %if
            tic; % starts a stopwatch timer
            [x(m+1:msp), error, iter, iflag, res_iter, m_v] = bicg2psx(A(:,m+1:msp), x(m+1:msp), b(m+1:msp), max_it, tol,precond);
            tempo=toc; operaz=flops;
        elseif (mspl==5) % CGN
            if (precond<=0)
                fprintf('\nunprecond. CG-normal eq.\n')
            else
                fprintf('\nLeft precond. CG-normal eq.\n')
            end %if
            tic; % starts a stopwatch timer
            [x(m+1:msp), error, iter, iflag, m_v]=cgne(A(:,m+1:msp), x(m+1:msp), b(m+1:msp), max_it, tol,precond);
            tempo=toc; operaz=flops;
        elseif (mspl==8)
            fprintf(1,'\n iterazioni con la sola matrice circolante \n')
            tic; % starts a stopwatch timer
            x(m+1:msp)=dbprecb(b(m+1:msp)); iter=0; iflag=0; res_iter=0; error=0;
            tempo=toc; operaz=flops; m_v=0;
        elseif (mspl==9) % direct factorization
            fprintf(1,'\nMatlab b.u. Gaussian Elimination ...\n')
            tic; % starts a stopwatch timer
            x(m+1:msp)=A(:,m+1:msp)\b(m+1:msp); % soluzione diretta
            tempo=toc; operaz=flops; iflag=0; iter=0; error=-1; res_iter=[]; m_v=0;
        elseif (mspl==10)
            tic; % starts a stopwatch timer
            [x(m+1:msp),iflag,error,iter,res_iter] = idrs(@(x) Mprod(A1,B1,L,hi,x),b(m+1:msp),4,tol,max_it,@(v) dbprec(v),[],x(m+1:msp));
            tempo=toc; operaz=flops; m_v = 3*iter;
        elseif (mspl==11) % FGMRES+dbprec(v) with variable solver with Recycling Subspace Method
            %fprintf(1,'\nBiCG...\n')
            tic; % starts a stopwatch timer
            [x(m+1:msp),iflag,res,iter] = fgmres(@(x) Mprod(A1,B1,L,hi,x),b(m+1:msp),1000,max_it,tol);
            tempo=toc; operaz=flops;
            if iflag ~= 0
                iter = -iter;
            end
            m_v = iter;
        end %if mspl

        x(1:m)=y0;

        if (0)
            if (mspl~=10 | mspl~=2)
                if iflag~=0
                    %disp('Grave problema nella soluzione del sist. lin.');
                end %if flag==0
                %fprintf('\ntime=%f Mflops=%i\n',tempo,operaz/1000000)
            end %if
        end %if
        if iflag~=0
            %disp('Problema nella soluzione del sist. lin.');
        end %if flag==0
        %fprintf('\ntime=%f Mflops=%i\n M-v=%i\n',tempo,operaz/1000000,m_v)
        %A=A(:,m+1:msp); % 'squadra' A togliendo la prima block-colonna
        fprintf('\n%i & %i & %i & %f',N,s,m_v,tempo);
        keycl=2;
        % keycl=input('[0] plot u(x,t) [1] restart [2] end:');
        if keycl==1
            itlmf
        elseif keycl==2
            % do nothing, esce
        else
            % FUNZIONI ACCESSORIE:
            % Grafico delle soluzioni, ecc.

            figure(2)

            tx(1)=t0;
            y1(1,:)=y0;
            for i=2:(s*p+1)
                %tx(i)=tx(i-1)+h_tilde(ceil((i-1)/s))/s;
                tx(i)=tx(i-1)+hi; % passo costante esterno
            end
            for i=2:(s*p+1)
                for c=1:m
                    y1(i,c)=x(m*(i-1)+c);
                end
            end

            if (pb>9 && pb~=13 && pb~=19 ) % Rimosso da Fabio
                if (pb~=11 && pb~=15)
                    xx=dx:dx:(x_max-dx); tx=0:hi:T;
                else
                    xx=0:dx:(x_max-dx); tx=0:hi:T;
                end %if
                %clf ***
                [Xm,Ym]=meshgrid(xx,tx);
                if pb == 21 || pb == 22 || pb == 23 || pb == 24 % Aggiunto da Fabio
                    subplot(1,2,1); surf(Xm,Ym,y1(:,1:N-1));
                else
                    surf(Xm,Ym,y1(:,1:N));
                end
                view(120,30) %view(20,30); % vista come \cite{libroBT}=view(120,30)
                %surfl(Xm,Ym,y1(:,1:N));
                %colormap(gray)
                %shading interp
                xlabel('x'); ylabel('t'); zlabel('u(x,t)');
                file=[titolo ' (h=', num2str(hi) ,')'];
                title(file,'Interpreter','LaTeX')
                if pb == 23
                    for j=1:length(tx)
                        ytrue(j,:) = 4*exp(-tx(j)).*(xx.^2).*(2-xx).^2;
                    end
                    subplot(1,2,2); surf(Xm,Ym,ytrue);
                    file=[titolo '(True Solution)'];
                    title(file,'Interpreter','LaTeX');
                    xlabel('x'); ylabel('t'); zlabel('u(x,t)');
                    view(120,30)
                    figure(3)
                    surf(Xm,Ym,abs(y1(:,1:N-1)-ytrue));
                    xlabel('x'); ylabel('t'); zlabel('|u_e(x,t)-u(x,t)|');
                    view(120,30)
                end
            elseif (pb<10 || pb==19 ) % Rimosso da Fabio: || pb==22

                clf
                subplot(2,1,1)
                hold on
                for i=1:m
                    if (m<9)
                        plot(tx,full( y1(:,i) ),color(i));
                    else
                        plot(tx,full( y1(:,i) ),'g');
                    end %if
                end %i
                quota=axis;
                db=ones(1,s*p+1)*quota(3);
                plot(tx,db,'+')
                file=['(a) t  (h=', num2str(hi) ,')'];
                xlabel(file)
                ylabel('Y(t)')
                quota=axis;
                axis([t0 T quota(3) quota(4)])
                hold off
                if (mspl~=8)
                    %figure(1)
                    subplot(2,1,2)
                    semilogy(1:length(res_iter),res_iter,'g')
                    % segue: numero prodotti M-v
                    % (n=numero M*v nell'algoritmo Bi-CG)
                    %semilogy(n:n:n*length(res_iter),res_iter,'g') % per n prodotti M-v
                    ylabel('||r||')
                    xlabel('Iterations')
                    %quota=axis;
                    %axis([1 iter error 1e+1])
                else
                    % do nothing, usare programma
                    % esterno per grafico
                end %if

                % GRAFICO ERRORI GLOBALI PER GALLERY


                subplot(3,1,3)
                if pb==2
                    for i=1:s*p+1  % per problema 2 (oscillatore armonico)
                        txs(1)=tx(i);
                        sqrtg=sqrt(gamm);
                        yvero(i,1)=cos(sqrtg*txs); % per problema n.2 (m=2)
                        yvero(i,2)=-sqrtg*sin(sqrtg*txs);
                    end
                    plot_error=1;
                elseif pb==4
                    for i=1:s*p+1  % per problema 4
                        txs(1)=tx(i);
                        yvero(i,1)=2*exp(-txs)+sin(txs);
                        yvero(i,2)=2*exp(-txs)+cos(txs);
                    end
                    plot_error=1;
                    %elseif pb==5
                    %  for i=1:s*p+1  % per problema 5
                    %     txs(1)=tx(i);
                    %    yvero(i,1)=exp(-txs)-exp(-300*txs);
                    %    yvero(i,2)=exp(-txs)+exp(-300*txs);
                    %  end
                    %  plot_error=1;
                elseif pb==1
                    % per problema 2 (oscillatore armonico)
                    %txs(1)=tx(i);
                    yvero=(exp(-tx)-exp(-1e+5*tx))';
                    plot_error=1;
                elseif pb==3
                    for i=1:s*p+1  % per problema 3 (mildly stiff)
                        txs(1)=tx(i);
                        yvero(i,1)=(exp(-2*txs)+exp(-40*txs)*( cos(40*txs)+sin(40*txs) ) )/2; % per problema 3 (m=3)
                        yvero(i,2)=(exp(-2*txs)-exp(-40*txs)*( cos(40*txs)+sin(40*txs) ) )/2;
                        yvero(i,3)=             exp(-40*txs)*(-cos(40*txs)+sin(40*txs));
                    end
                    plot_error=1;
                elseif pb==5
                    %for i=1:s*p+1  % per problema 5 (instabile)
                    %  txs(1)=tx(i);
                    yvero=(sin(tx)+cos(tx))';
                    %end
                    plot_error=1;
                elseif pb==7
                    for i=1:s*p+1  % per problema 3 (mildly stiff)
                        txs(1)=tx(i);
                        yvero(i,1)=exp(-50*txs)+exp(0.1*txs)*sin(8*txs); % per problema 7 (m=3)
                        yvero(i,2)=-exp(-50*txs)+exp(0.1*txs)*cos(8*txs);
                        yvero(i,3)=exp(-50*txs)+exp(0.1*txs)*(cos(8*txs)+sin(8*txs));
                    end
                    plot_error=1;
                else
                    % do nothing
                end

                %plot_error=0;
                if (plot_error)
                    errore=yvero-y1;
                    if (mspl~=8 && pb<10)
                        figure(1)
                    else
                        subplot(2,1,2)
                    end %if
                    semilogy(tx,abs(errore));
                    xlabel('t  Global error')
                end % if

            end %if pb>9

        end %if keycl

        % fine itlmf.m

        clearvars -except Nvalue svalue precond pb pol_flag polflag metodo ...
            passi restart_it fracprec diagNdiv fracval apprJ mspl
    end
%end
fprintf('\n');