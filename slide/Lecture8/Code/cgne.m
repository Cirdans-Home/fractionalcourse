function x = cgne ( A, b, x )

%*****************************************************************************80
%
%% cgne() carries out the conjugate gradient method for the normal equations.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 March 2020
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    real A(n,n): the matrix.
%
%    real B(n,1): the right hand side.
%
%    real X(n,1): an initial guess for the solution.
%
%  Output:
%
%    real X(n,1): an improved solution estimate.
%
  n = size ( A, 1 );

  r = b - A * x;
  z = A' * r;
  d = z;

  for i = 1 : n

    alpha = ( z' * z ) / ( ( A * d )' * ( A * d ) );

    x = x + alpha * d;

    r = b - A * x;
    if ( norm ( r ) == 0.0 )
      fprintf ( 1, '  Early finish on step %d\n', i );
      break
    end

    zold = z;
    z = A' * r;

    beta = ( z' * z ) / ( zold' * zold );

    d = z + beta * d;

  end

  return
end

