function I = constfracint(f,a,t,alpha,N,lambda)
%CONSTFRACINT computes the fractional integral with the weighted piecewise
%constant approximation of the function f between a and t, over N uniformly
%distributed intervals.
h  = (t-a)/N;
tk = (a:h:t)';
b  = zeros(N,1);
for k=0:N-1
   b(k+1) = (k+1)^alpha - k^alpha; 
end
b = h^alpha*b/gamma(alpha+1);
p = f(tk);
I = flipud(b)'*(lambda*p(1:N) + (1-lambda)*p(2:N+1));
end

