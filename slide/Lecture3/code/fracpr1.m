function [t,y] = fracpr1(f,tspan,y0,alpha,tau)
%FRACPR1 Uses the Fractional Product Integral Rule based on the rectangule
%formula to solve for: cD_{tspan}^\alpha y(t) = f(t,y(t)) with y(0) = y0
%   INPUT: f = @(t,y) f(t,y) 
%          tspan = [t0,tf] 
%          y0 = initial condiition
%          alpha = order of the fractional derivative
%          tau =  time-step (uniform)

t = tspan(1):tau:tspan(2);
N = length(t);
b = zeros(N,1);
if isrow(y0)
   y0 = y0';
end
k = length(y0);
y  = zeros(k,N);
fy = zeros(k,N);

for i=0:N-1
% These are the values b0,b1,...,bN-1 we need for all the computations
   b(i+1) = ((i+1)^alpha - i^alpha)/alpha; 
end
y(:,1) = y0;
c = (tau^alpha/gamma(alpha));
for j=1:N-1
    fy(:,j) = f(t(j),y(:,j)); % New value of f
    y(:,j+1) = y(:,1) + c*sum(flipud(b(1:j)).'.*fy(:,1:j));
end

end

