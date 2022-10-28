function [ro,si] = coef( type, k, i )
%
%  Coefficienti di TOM / ETR (GAM) / ETR2 (GAM2) /GBDF.
%
if nargin==2    % TOM, che non richiede il parametro I, TYPE non utilizzato;
  [ro,si] = rosi_tom( k );
elseif type==1  % ETR del primo tipo, o GAM;
  [ro,si] = rosi_gam( k, i );
elseif type==2  % ETR del secondo tipo, o GAM2;
  [ro,si] = rosi_etr( k, i );
elseif type==0  % BDF generalizzate.
  [ro,si] = rosi_bdf( k, i );
elseif type==-1  % midpoint
  [ro,si] = rosi_bdf( k, i );
else
  [ro,si] = rosi_et2( k, i );
end
