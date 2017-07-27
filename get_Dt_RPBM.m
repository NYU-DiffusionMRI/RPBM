function Dt=get_Dt_RPBM(t,zeta)
% Finds D(t/tau) from D(omega) from Nature Physics 7:508 (2011), Eq (6)
% by deforming the contour of integration from the Re omega axis 
% to the one around the negative imaginary axis in the lower half plane.
% A few terms are subtracted so that the integrand is regular at zero.
%
% (c) Dmitry S Novikov 2015 dima@alum.mit.edu

lambda=1+zeta; % tortuosity
A=2*(sqrt(lambda)-1)/lambda^2;
B=2*(sqrt(lambda)-1)*(sqrt(lambda)-2)/lambda^3;
C=(8*(sqrt(lambda)-1)^2-zeta*sqrt(lambda))/lambda^4;

fint2=inline('(A*sqrt(y)+C*y.^(3/2)+imag(1./(1+zeta+2*i*sqrt(y).*(i*sqrt(y)-1).*(1 - sqrt(1 + zeta./(i*sqrt(y)-1).^2))))).*exp(-t*y)./y.^2', 'y', 't', 'zeta', 'A', 'C');

%ymin=1e-10; ymax=1e10; tol=1e-10;
Dt=integral(@(y)fint2(y,t,zeta,A,C), 0, Inf, 'ArrayValued', true)./(pi*t) + 1/lambda + 2*A./sqrt(pi*t) + B./t - C*t.^(-3/2)/sqrt(pi);
