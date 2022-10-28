function v = Psolve(omega,phi,psi,h,L,v)
%PSOLVE Apply the structured preconditioner

m = length(psi);    % Number of time steps
n = size(L,1);      % Number of space dofs
Omega = spdiags(omega.^(-(0:m-1).'./m),0,m,m);
I = speye(n,n);

v = reshape(v,n,m);
% Step 1
v = -fft(v/Omega',[],2); % F^T = F
% Step 2: The should be parallel part
for i=1:m
    v(:,i) = (phi(i)*I - h*psi(i)*L)\v(:,i);
end
% Step 3
v = -ifft(v,[],2)/Omega;
v = real(v(:));

end