function [jacobian] = getJacobianNonDimensionalParFor2(X,p)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [jacobian] = getJacobianNonDimensionalParFor2(X,p)
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
% Note: this file can be easily adapted to use Matlab's parallelization 
% toolbox as indicated in the loop below.
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

dQ2dq2 = ep;
dQ2dc = -1;
dQ1dc = 1./Q1x.*(rho*Q2x.*dQ2dc - (rho-1)*c*(1 + ep^2*mu^2*etaX.^2));
dQ2dEta = -w2*ep;
dQ1dq2 = 1./Q1x.*(rho*Q2x)*ep;


jacobian = zeros(2*M+1,2*M+1);
jacobianA = zeros(M,M); jacobianB = zeros(M,M);
jacobianC = zeros(M,M); jacobianD = zeros(M,M);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Change the "for" to "parfor" to use Matlab's parallelization toobox.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for iter = 1:M^2
    nInd = mod(iter,M);
    if nInd==0
        nInd = M;
    end
    
    n = nInd - N - 1; k = n;
    
    mInd = (iter - nInd)/M+1; m = mInd - N - 1; l = m;
    
    coshTerm1 = cosh(mu*k*ep*eta) - tanh(mu*k)*sinh(mu*k*ep*eta);
    sinhTerm1 = sinh(mu*k*ep*eta) - tanh(mu*k)*cosh(mu*k*ep*eta);
    coshTerm2 = cosh(mu*k*ep*eta) + tanh(mu*k*del)*sinh(mu*k*ep*eta);
    sinhTerm2 = sinh(mu*k*ep*eta) + tanh(mu*k*del)*cosh(mu*k*ep*eta);
    
    if mInd==1
        jacobianE(iter,1) = hat(sinhTerm1.*dQ1dc,n);%
        jacobianF(iter,1) = hat(sinhTerm2.*dQ2dc,n);%
    end
    
    dQ1dEta = 1./Q1x.*(rho*Q2x*dQ2dEta - (rho-1)*(1i*l*etaX.*(c^2-2*ep*eta)*ep^2*mu^2 - ep*(1 + ep^2*mu^2*etaX.^2)));
    dAdEta = w1*sinhTerm1*ep + (mu*k)*coshTerm1.*Q1x*ep + sinhTerm1.*dQ1dEta;
    dBdEta = w2*sinhTerm2*ep + (mu*k)*coshTerm2.*Q2x*ep + sinhTerm2.*dQ2dEta;
    dAdq = sinhTerm1.*dQ1dq2;
    dBdq = sinhTerm2.*dQ2dq2;
    
    if abs(n-m)<=N
        jacobianA(iter) = (abs(n)>0)*(abs(m)>0)*hat(dAdEta,n-m);%
        jacobianB(iter) = (abs(n)>0)*(abs(m)>0)*hat(dAdq,n-m);%
        jacobianC(iter) = (abs(n)>0)*(abs(m)>0)*hat(dBdEta,n-m);%
        jacobianD(iter) = (abs(n)>0)*(abs(m)>0)*hat(dBdq,n-m);%
    end
    
end

jacobianA = reshape(jacobianA,[M M]);
jacobianB = reshape(jacobianB,[M M]);
jacobianC = reshape(jacobianC,[M M]);
jacobianD = reshape(jacobianD,[M M]);
jacobian(1:2*M,1:(2*M+1)) = [jacobianA jacobianB jacobianE; jacobianC jacobianD jacobianF];

jacobian(N+1,N+1) = 1;
jacobian(M+N+1,M+N+1) = 1;


switch p.ContinuationType% 1 = FourierMode, 2 = InfNorm, 3 = Max Value
    case 1
        jacobian(end,[N N+2]) = 1;
    case 2
        jacobian(end,1:M) = (ones(1,M));
    case 3
        jacobian(end,1:M) = .5*(ones(1,M) - (-1).^((-N:N)));
        
end
