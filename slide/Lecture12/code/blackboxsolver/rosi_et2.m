function [ro,si] = rosi_et2( k, j )
%
% function [ro,si] = rosi_et2( k, j )
%
% Coefficienti dell' ETR2s (generalizzato) a K passi e J condizioni iniziali.
% Quando k=2 nu-1, e j=nu, si ottengono gli ETR2s veri e propri.
% Polinomio sigma fissato.
%
ro = zeros( k+1, 1 ); si = ro; s  = 1; 
for i = 2:k+1; 
  ro(i) = s*(i-1); s = -s; 
end
ro = .5*ro; ro(2) = 1;
ro = vsolve( -j:k-j, ro );
si(j:j+1) = [.5; .5];
