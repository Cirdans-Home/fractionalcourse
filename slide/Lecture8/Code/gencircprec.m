%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : gencircprec.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Circulant Preconditioner Apply O(n log(n))
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ c ] = gencircprec( t,type )
%GENCIRCPREC This function generates the circulant preconditioner relative
%to the Hermitian Toeplitz matrix with first column t.
%   TYPE OF PRECONDITIONER:
%   'tchan'     Tony Chan Preconditioner
%   'strang'    Strang Preconditioner
%   'rchan'     Raymond Chan Preconditioner
%   'dirichlet' Preconditioner from modified Dirichlet kernel
%   'poussin'   Preconditioner from de la Vallee Poussin kernel
%   'hann'      Preconditioner from von Hann kernel
%   'hamming'   Preconditioner from Hamming kernel
%   'bernstein' Preconditioner from Bernstein kernel
%   'superoptimal' Superoptimal preconditioner

n = size(t,1);          % size of the preconditioner
t1 = conj(t(n:-1:2));   % is the last column of T_n(f)

if strcmp(type,'tchan')         % T. Chan’s preconditioner
    coef = (1/n:1/n:1-1/n)';
    c = fft([t(1);(1-coef).*t(2:n)+coef.*t1]);
elseif strcmp(type,'strang')    % Strang’s preconditioner.
    c = fft([t(1:n/2);0;conj(t(n/2:-1:2))].')';
elseif strcmp(type,'rchan')     % R. Chan’s preconditioner
    c = fft([t(1);t(2:n)+t1].')';
elseif strcmp(type,'dirichlet') % Modified Dirichlet kernel
    coef = [t(1);t(2:n)+t1].';
    coef(2) = coef(2)-0.5*t1(1);
    coef(n) = coef(n)-0.5*t(n);
    c = fft(coef)';
elseif strcmp(type,'poussin')   % de la Vallee Poussin kernel
    coef = zeros(1,n);
    coef(1) = t(1);
    m = floor(n/2);
    coef(2:m+1) = (1:m).*conj(t(2*m:-1:m+1)')/m+t(2:m+1)';
    coef(m+2:2*m) = (m-1:-1:1).*t(m+2:2*m)'/m+conj(t(m:-1:2)');
    c = fft(coef)';
elseif strcmp(type,'hann')      % von Hann kernel
    tp = pi/(2*n);
    coef = ((cos(tp:tp:(n-1)*tp)).^2)';
    coef = [t(1);coef.*t(2:n)+(1-coef).*t1];
    c = fft(coef);
elseif strcmp(type,'hamming')   % Hamming kernel
    tp = pi/(2*n);
    coef = ((cos(tp:tp:(n-1)*tp)).^2)';
    c = [t(1);0.46*(coef.*t(2:n)+(1-coef).*t1)+0.54*(t(2:n)+t1)];
    c = fft(c);
elseif strcmp(type,'bernstein') % Bernstein kernel
    tp = pi/n;
    coef = (0.5*(1+exp((tp:tp:(n-1)*tp)*1i)))';
    c = [t(1);coef.*t(2:n)+(1-coef).*t1];
    c = fft(c);
elseif strcmp(type,'superoptimal')  % Superoptimal preconditioner
    h = zeros(n,1); % h: first column of the circulant part
    s = zeros(n,1); % s: first column of the skew-circulant part
    h(1) = 0.5*t(1);
    s(1) = h(1);
    h(2:n) = 0.5*(t(2:n) + t1);
    s(2:n) = t(2:n)-h(2:n);
    ev1 = fft(h);
    coef = (1:-2/n:2/n-1)';
    c = coef.*s;                % first column of T. Chan’s preconditioner
                                % for the skew-circulant part
    ev2 = fft(c);
    d = (exp((0:1/n:1-1/n)*pi*1i)).';
    s = s.*d;
    sev = fft(s);               % eigenvalues of the skew-circulant part
    sev = sev.*conj(sev);
    s = ifft(sev);
    s = conj(d).*s;
    h = coef.*s;
    ev3 = fft(h);
    c = (abs(ev1).^2 + 2*ev1.*ev2 + ev3)./ev1;
else
    fprintf('Preconditioner not recognized, give back identity\n');
    c = ones(n,1);
end

c = real(c);

end

