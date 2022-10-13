function [g] = gl(n,alpha)
%GL Produces the N+1 Grunwald-Letnikov coefficients for a given alpha
g = cumprod([1, 1 - ((alpha+1) ./ (1:n))]);
end

