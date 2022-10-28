function [a,b] = mab( type, k, n )
%
% type: 0 (GBDF)
%       1 (GAM)
%       2 (ETR2)
%       3 (TOM)  
%       4 (ETR2/1)
%
%  ATTENZIONE: E' STATA RIMOSSA LA RIGA CHE AGGIUNGE LA MATRICE IDENTITA' AD A ECC.
%
if type==3
  if nargin==2, n = 2*k-1; end
  nu = ceil( k/2 ); 
  if k-2*nu+1>0
    disp('SOLO NUMERO DI PASSI DISPARI')
    return, 
  end  % k dispari solamente
  
  a = spalloc( n, n+1, (k+1)*n+k^2 ); b = a;
  for i = 1:nu-1
    [ro,si] = coef( 1, 2*k-1, i );
    a(i,1:2*k) = ro.';
    b(i,1:2*k) = si.';
  end
  [ro,si] = coef( 3, k ); ro = ro.'; si = si.';
  for i = nu:n-nu+1
    a(i,i-nu+(1:k+1)) = ro;
    b(i,i-nu+(1:k+1)) = si;
  end
  ini = 2*k-1;
  for i = n:-1:n-nu+2
    [ro,si] = coef( 1, 2*k-1, ini );
    a(i,n+1-2*k+(1:2*k)) = ro.';
    b(i,n+1-2*k+(1:2*k)) = si.';
    ini = ini-1;
  end
else
  if nargin==2, n = k+(type>3); end
  if type==0
    nu = fix( (k+2)/2 );
  elseif type==-1
    nu = ceil( k/2 );
  else
    nu = ceil( k/2 );
  end
  a  = spalloc( n, n+1, (k+1)*n ); b = a;
  for i = 1:nu
    plus1 = type>3 & i<nu;
    [ro,si] = coef( type, k+plus1, i );
    a(i,1:k+1+plus1) = ro.';
    b(i,1:k+1+plus1) = si.';
  end
  for i = nu+1:n-(k-nu)
    a(i,i+1+(-nu:k-nu)) = ro.';
    b(i,i+1+(-nu:k-nu)) = si.';
  end
  plus1 = (type>3) ; j = nu + plus1;  
  for i = n-(k-nu)+1:n
    j = j + 1;
    [ro,si] = coef( type, k+plus1, j ); 
    a(i,n+1+(-k-plus1:0)) = ro.';
    b(i,n+1+(-k-plus1:0)) = si.';
  end
end
%[m,n] = size(a);
% a = [1 zeros(1,m); a]; b = [zeros(1,n); b];