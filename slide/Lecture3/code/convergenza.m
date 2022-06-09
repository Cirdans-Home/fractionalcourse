function q = convergenza(e)
%%CONVERGENZA produce una stima dell'ordine di convergenza della
%%successione x_n ad xtrue.
  q = zeros(length(e)-1,1);
  for n = 2:(length(e))
      q(n-1) = log(e(n+1)/e(n))/log(e(n)/e(n-1));       
  end
end