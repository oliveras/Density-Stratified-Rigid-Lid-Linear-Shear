function [L11 L11m L12m L13m R11om R11 R12 R13] = getConstantMatricesNonDimensional(X, p)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [L11 L11m L12m L13m R11om R11 R12 R13] = getConstantMatricesNonDimensional(X, p)
%
% Computes the spectral matrices that are multiplied directly by the 
% Floquet parameter (FP) in order to simplify and speed up the computation 
% of the spectra.  Note: these terms need only to be computed once for all
% Floquet parameters. All of these terms come from the linearization of the
% Bernoulli (or pressure) equation.
%
% INPUT
%     - 'X' consists of the Fourier coefficients of eta, followed by the
%       Fourier coefficients of q2x, and then the wave speed c.
%     - 'p' consist of model parameters passed through as a struct.
%     
% OUTPUT
%      - L11 - The component of L11 not multiplied by the FP
%      - L11m - The component of L11 multiplied by (dx + 1i*FP)
%      - L12m - The component of L12 multiplied by (dx + 1i*FP)
%      - L13m - The component of L13 multiplied by (dx + 1i*FP)
%      - R11om - The component of R11 'divided' by (dx + 1i*FP)
%      - R11 - The component of R11 not multiplied by the FP
%      - R12 - The component of R12 not multiplied by the FP 
%      - R13 - The component of R13 not multiplied by the FP
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = p.N;                    M = p.M;                L = p.L;
x = (0:2*N)'/M*L;           dx = x(2)-x(1);         qSign=p.qSign;
w1 = p.w1;                  w2 = p.w2;              rho = p.rho;
ep = p.epsilon;             mu = p.mu;              del = p.delta;
numModes = p.numModes;
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

Omega = w1 - rho*w2;
T1 = (Q1x)./(1+(ep*mu*etaX).^2);
T2 = rho*(Q2x)./(1+(ep*mu*etaX).^2);


zM = zeros(2*numModes+1,2*numModes+1);

L11 = zM; L11m = zM;
L12m = zM; L13m = zM;
R12 = zM; R13 = zM;
R11om = zM; R11 = zM;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R1
% - - - - - - - - - -
% R11 dx^0 term
R11h = ep*mu^2*hat((T2 - T1).*etaX);
% R11 dx^(-1) term
R11om = -Omega*eye(size(zM));                        
% - - - - - - - - - -
R12 = eye(size(zM)); 
% - - - - - - - - - -
R13 = -rho*eye(size(zM));           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L11h = hat(rho - 1 + T1*w1 - T2*w2);
L11mh = ep*mu^2*hat((T1.^2 - T2.^2/rho).*etaX);
L12mh = -hat(T1);                             % STerm dx^1
L13mh = hat(T2);                            % TTerm dx^2

for m = -numModes:numModes
    for n = -numModes:numModes
        mInd = numModes + m + 1;
        nInd = numModes + n + 1;
        if abs(m-n)<=numModes
            L11(mInd,nInd) = L11h(m-n+N+1);
            L11m(mInd,nInd) = L11mh(m-n+N+1);
            L12m(mInd,nInd) = L12mh(m-n+N+1);
            L13m(mInd,nInd) = L13mh(m-n+N+1);
            R11(mInd,nInd) = R11h(m-n+N+1);
            
        end
    end
end

