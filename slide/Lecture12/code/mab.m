function [a,b] = mab( k, n )
nu = fix( (k+2)/2 );
a  = spalloc( n, n+1, (k+1)*n ); b = a;
for i = 1:nu
    [ro,si] = rosi_bdf( k, i );
    a(i,1:k+1) = ro.';
    b(i,1:k+1) = si.';
end
for i = nu+1:n-(k-nu)
    a(i,i+1+(-nu:k-nu)) = ro.';
    b(i,i+1+(-nu:k-nu)) = si.';
end
j = nu;
for i = n-(k-nu)+1:n
    j = j + 1;
    [ro,si] = rosi_bdf( k, j );
    a(i,n+1+(-k:0)) = ro.';
    b(i,n+1+(-k:0)) = si.';
end
end