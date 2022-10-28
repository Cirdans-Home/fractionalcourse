function [ro,si] = rosi_bdf( k, j )
%
%  [RO,SI] = ROSI_BDF( K, J )    Calcola i polinomi RO e SIGMA della BDF
%                                generalizzata (GBDF) a K passi e J con-
%                                dizioni iniziali (J=K sono le BDF IVM,
%                                J=fix(K/2)+1 da le GBDF A-stabili).
%  Coefficienti dei polinomi per potenze crescenti.
%
b  = zeros(k+1,1); b(2) = 1;
ro = vsolve( -j:k-j, b(:) );
si = zeros( k+1, 1 ); si( j+1 ) = 1;
