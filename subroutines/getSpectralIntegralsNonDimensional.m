function [intOut jacobian] = getSpectralIntegralsNonDimensionalN(X,p,options)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [intOut jacobian] = getSpectralIntegralsNonDimensionalN(X,p,options)
%
%
%
% Calculates the integral relations in order to solve for traveling wave
% solutions in the form F(X) = 0.
%
% INPUT
%     - 'X' consists of the Fourier coefficients of eta, followed by the
%       Fourier coefficients of q2x, and then the wave speed c.
%     - 'p' consist of model parameters passed through as a struct.
%     - 'options' options for using Matlab's fsolve command (part of the
%       optimization toolbox).
%
% OUTPUT
%      - 'intOut' consists of the left-hand side of the equation F(X) = 0
%      - 'jacobian' is the Jacobian for use in the nonlinear solver.
%
% Note: the Jacobian is calculated using the external file
% 'getJacobianNonDimensionalParFor2' which can be easily adapted to use 
% Matlab's parallelization toolbox as indicated in the file.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin==2
    options = optimset('Jacobian','off');
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = p.N;                    M = p.M;                L = p.L;
x = (0:2*N)'/M*L;           dx = x(2)-x(1);         qSign=p.qSign;
w1 = p.w1;                  w2 = p.w2;              rho = p.rho;
ep = p.epsilon;             mu = p.mu;              del = p.delta;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
etaHat = X(1:M); q2xHat = X((M+1):(2*M)); c = (X(end));

eta = invHat(etaHat);       etaX = d(eta);          q2x = invHat(q2xHat);
eta = real(eta);            etaX = real(etaX);      q2x = real(q2x);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Note: QjX = ep*qjx - wj*ep*eta - c;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Q2x = ep*q2x - w2*ep*eta - c;
Q1x = (qSign*sqrt(rho*(Q2x).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + (ep*mu*etaX).^2)));

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


for nInd = 1:M
    n = nInd - N - 1;
    k = n;
    
    coshTerm1 = cosh(mu*k*ep*eta) - tanh(mu*k)*sinh(mu*k*ep*eta);
    sinhTerm1 = sinh(mu*k*ep*eta) - tanh(mu*k)*cosh(mu*k*ep*eta);
    coshTerm2 = cosh(mu*k*ep*eta) + tanh(mu*k*del)*sinh(mu*k*ep*eta);
    sinhTerm2 = sinh(mu*k*ep*eta) + tanh(mu*k*del)*cosh(mu*k*ep*eta);
    
    integrand1 =  w1 * coshTerm1 / (mu * k)  + sinhTerm1 .* Q1x ;
    integrand2 =  w2 * coshTerm2 / (mu * k)  + sinhTerm2 .* Q2x ;
    
    intOut1(nInd,1) = hat(integrand1,n);
    intOut2(nInd,1) = hat(integrand2,n);
    
    if n==0
        intOut1(nInd,1) = etaHat(nInd);
        intOut2(nInd,1) = q2xHat(nInd);
    end
end

switch p.ContinuationType % 1 = FourierMode, 2 = InfNorm, 3 = Max Value
    case 1
        ampCond = etaHat(N) + etaHat(N+1) - p.amp/2;
    case 2
        ampCond = max((eta))-p.amp;
    case 3
        ampCond = 1/2*sum((ones(M,1) - (-1).^((-N:N)')).*etaHat)-p.amp;
end

intOut = real([intOut1;intOut2; ampCond]);

if strcmp(options.Jacobian,'off')
    jacobian=[];
else
    jacobian = getJacobianNonDimensionalParFor2(X,p);
    jacobian= real(jacobian);
    
end

jacobian = (jacobian);
intOut = (intOut);
