function f = vsolve( x, b )
%
%     f = vsolve( x, b )      Risolve il sistema lineare W(x) f = b,
%                             dove  W(x)  e' la Vandermonde definita
%                             dagli elementi del vettore x.
%
f = b;
n = length( x )-1;
for k = 1:n
  for i = n+1:-1:k+1
    f(i) = f(i) - x(k)*f(i-1);
  end
end
for k = n:-1:1
  for i = k+1:n+1
    f(i) = f(i)/( x(i) - x(i-k) );
  end
  for i = k:n
    f(i) = f(i) - f(i+1);
  end
end
