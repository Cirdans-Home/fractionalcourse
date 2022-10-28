function [ro,si] = rosi_etr( k, j )
%
% Coefficienti dell' ETR2s (generalizzato) a K passi
% e J condizioni iniziali.
% Quando k=2 nu-1, e j=nu, si ottengono gli ETR2s veri e propri.
% Coefficienti dei polinomi per potenze crescenti.
%
a1 = zeros(k+1,1); a2 = a1; si = a1;
a1(2) = 1; a2(3:k+1) = [2:k]'.*cumprod( -ones(k-1,1) );
vt = [-j:k-j].^(k+1);
a1 = vsolve( -j:k-j, a1 );
a2 = vsolve( -j:k-j, a2 );
beta = ( (k+1)*(-1)^k -vt*(a2+a1) )/( (k+1)*(-1)^k -vt*a2 );
ro   = a1 + (1-beta)*a2;
si(j) = 1-beta; si(j+1) = beta;
