function [Psi, DPsi] = monomial_basis_sin(deg, dim)
%  [Psi, DPsi] = monomial_basis(deg, dim) returns a monomial basis function and its derivative
%   There is not a 1 included in Psi, only the linear and higher order terms
%   deg = degree of monomial
%   dim = number of states
k = linspace(2, deg, deg-1);
d = dim;
x=sym('x',[d,1]);
assume(x,'real')

Psi = [x.'];
for i=1:size(k,2)
    m = nchoosek(k(i)+d-1,d-1);
    dividers = [zeros(m,1),nchoosek((1:(k(i)+d-1))',d-1),ones(m,1)*(k(i)+d)];
    a = diff(dividers,1,2)-1;
    for i = 1:size(a,1)
        Psi = [Psi prod(x.' .^ a(i,:))];
    end
end
DPsi = jacobian(Psi,x);
Psi = Psi'; % monomial only

x=sym('x',[d*2,1]);

Psi = [x(2); Psi]
end
%   Psi = [Psi;sin(x(1));cos(x(1)); x(1).*sin(x(1));x(2).*sin(x(1));x(1).*cos(x(1));x(2).*cos(x(1))];
%     Psi = [Psi; x(1).*sin(x(1));x(2).*sin(x(1));x(1).*cos(x(1));x(2).*cos(x(1))];
%     Psi = [Psi; sin(x(1))-x(1); x(2).*sin(x(1));(x(2).*sin(x(1)))^2; x(2).*cos(x(1))-x(2); (x(2).*cos(x(1))-x(2))^2];
