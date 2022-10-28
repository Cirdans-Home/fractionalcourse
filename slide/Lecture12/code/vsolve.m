function f = vsolve( x, b )
%VSOLVE Vandermonde system solution
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
end