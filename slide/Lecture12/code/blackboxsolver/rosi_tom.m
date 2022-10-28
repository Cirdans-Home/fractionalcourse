function [ro,si] = rosi_tom( k )
%
%  [ro,si] = ROSI_TOM( K )	Coefficienti del TOM a K passi, K = 2*NU-1.
%                               Il metodo va usato con NU condizioni iniziali.
%  Coefficienti dei polinomi per potenze crescenti.
%
nu = ceil( k/2 );
a  = cvsolve( -nu:k-nu, [zeros(2*k+1,1); 1] );
si = -a(k+2:2*k+2); ro = a(1:k+1); ro = ro/sum(si); si = si/sum(si);
