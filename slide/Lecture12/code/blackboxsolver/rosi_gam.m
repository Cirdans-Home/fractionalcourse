function [ro,si] = rosi_gam( k, j );
%
%     [ro,si] = rosi_gam( k, j )       Polinomi RO e SI del metodo di
%                                     Adams generalizzato a K passi e
%                                     J condizioni iniziali. I metodi
%                                     di Adams-Moulton hanno J=K.
%                          Richiamato da GAMBLK.M; richiama VSOLVE.M.
%     Coefficienti dei polinomi per potenze crescenti.
%
if j<1 | j>k, error( 'J<K!' ), end
b  = zeros(k+1,1); s = 1; b(1) = 1;
for i = 2:k+1, s=-s; b(i) = s/i; end
si = vsolve( -j:k-j, b(:) );
ro = zeros( k+1, 1 ); ro(j:j+1) = [-1; 1];
if nargout==1, ro=si; end
