function [F,esterr] = funm2(A,fun)
%FUNM Evaluate general matrix function.
%   F = FUNM(A,FUN) for a square matrix argument A, evaluates the
%   matrix version of the function FUN. For matrix exponentials,
%   logarithms and square roots, use EXPM(A), LOGM(A) and SQRTM(A)
%   instead.
%
%   FUNM uses a potentially unstable algorithm.  If A is close to a
%   matrix with multiple eigenvalues and poorly conditioned eigenvectors,
%   FUNM may produce inaccurate results.  An attempt is made to detect
%   this situation and print a warning message.  The error detector is
%   sometimes too sensitive and a message is printed even though the
%   the computed result is accurate.
%
%   [F,ESTERR] = FUNM(A,FUN) does not print any message, but returns
%   a very rough estimate of the relative error in the computed result.
%
%   If A is symmetric or Hermitian, then its Schur form is diagonal and
%   FUNM is able to produce an accurate result.
%
%   L = LOGM(A) uses FUNM to do its computations, but it can get more
%   reliable error estimates by comparing EXPM(L) with A.
%   S = SQRTM(A) and E = EXPM(A) use completely different algorithms.
%
%   Example
%      FUN can be specified using @:
%         F = funm(magic(3),@sin)
%      is the matrix sine of the 3-by-3 magic matrix.
%
%   See also EXPM, SQRTM, LOGM, @.

%   C.B. Moler 12-2-85, 7-21-86, 7-11-92, 5-2-95.
%   Copyright 1984-2001 The MathWorks, Inc.
%   $Revision: 5.15 $  $Date: 2001/04/15 12:01:37 $

% Parlett's method.  See Golub and VanLoan (1983), p. 384.

if isstr(A), error('The first argument must be a matrix.'); end

[Q,T] = schur(A);
[Q,T] = rsf2csf(Q,T);
F = diag(feval(fun,diag(T)));
[n,n] = size(A);
dmin = abs(T(1,1));
for p = 1:n-1
   for i = 1:n-p
      j = i+p;
      s = T(i,j)*(F(j,j)-F(i,i));
      if p > 1
         k = i+1:j-1;
         s = s + T(i,k)*F(k,j) - F(i,k)*T(k,j);
      end
      d = T(j,j) - T(i,i);
      if d ~= 0 
         s = s/d;
      end
      F(i,j) = s;
      dmin = min(dmin,abs(d));
   end
end
F = Q*F*Q';
if isreal(A) && norm(imag(F),1) <= 10*n*eps*norm(F,1)
   F = real(F);
end

esterr=0;
%if dmin == 0, dmin = eps; end
%esterr = min(1,max(eps,(eps/dmin)*norm(triu(T,1),1)));
%if any(~isfinite(F)), esterr = Inf; end
%if (nargout < 2) & (esterr > 1000*eps)
%   disp(' ')
%   disp(['WARNING: Result from FUNM may be inaccurate.' ...
%          ' esterr = ' num2str(esterr)])
%end
