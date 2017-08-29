function [eVals ] = getSpectralMatricesNonDimensionalParFor(X,p,floquetParams)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [eVals ] = getSpectralMatricesNonDimensionalParFor(X,p,floquetParams)
%
% Calculates the spectrum associated with the traveling wave soltuion
% represented in vector X with associated parameters provided in p.
%
% INPUT
%     - 'X' consists of the Fourier coefficients of eta, followed by the
%       Fourier coefficients of q2x, and then the wave speed c.
%     - 'p' consist of model parameters passed through as a struct.
%     - 'floquetParams' consist of the Floquet parameters for which the
%        spectrum will be calculated.
%
% OUTPUT
%      - 'eVals' a matrix of eigenvalues where each column represents the
%         eigenvalues associated with a particular Floquet Parameter.  For
%         example,
%               eVals(j,:)
%         returns the eigenvalues associated with the j-th Floquet
%         parameter.
%
%      - 'eVecs' is the multi-dimensional matrix consisting of the
%         eigenvectors associated with particular Floquet parameters.  For
%         example,
%               squeeze(eVecs(j,:,:))
%         returns a matrix whose columns are the eigenvectors associated
%         with the j-th Floquet Parameter.
%
% Note: this file can be easily adapted to use Matlab's parallelization
% toolbox as indicated in the loop below.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = p.N;                    M = p.M;                L = p.L;
x = (0:2*N)'/M*L;           dx = x(2)-x(1);         qSign=p.qSign;
w1 = p.w1;                  w2 = p.w2;              rho = p.rho;
ep = p.epsilon;             mu = p.mu;              del = p.delta;
numModes = p.numModes;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% initialize the size of the eVals matrix.
eVals = zeros(2*(2*numModes+1),length(floquetParams));

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


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Compute the matrices that do not vary in the Floquet parameter
% separately.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[L11 L11m L12m L13m R11om R11 R12 R13] = getConstantMatricesNonDimensional(X, p);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Modify "for" to "parfor" to use the parallelization toolbox
for fInd = 1:length(floquetParams)
    RR21 = zeros(size(L11));
    zm = zeros(size(L11));
    LL21 = zm;
    LL22 = zm;
    RR31 = zm;
    LL31 = zm;
    LL33 = zm;
    f = floquetParams(fInd); DX = 1i*(f + meshgrid(-numModes:numModes));
    
    for m = -numModes:numModes
        k = f + m;
        
        coshTerm1 = (cosh(mu*ep*k*eta) - tanh(mu*k)*sinh(mu*ep*k*eta));
        sinhTerm1 = (sinh(mu*ep*k*eta) - tanh(k*mu)*cosh(mu*ep*k*eta));
        coshTerm2 = (cosh(mu*ep*k*eta) + tanh(k*mu*del)*sinh(mu*ep*k*eta));
        sinhTerm2 = (sinh(mu*ep*k*eta) + tanh(k*mu*del)*cosh(mu*ep*k*eta));
        
        
        
        
        coshTerm1Hat = hat(coshTerm1); sinhTerm1Hat = hat(sinhTerm1);
        coshTerm2Hat = hat(coshTerm2); sinhTerm2Hat = hat(sinhTerm2);
        
        coshTerm1Hat2 = hat(coshTerm1.*(Q1x));
        coshTerm2Hat2 = hat(coshTerm2.*(Q2x));
        
        for n = -numModes:numModes
            mInd = numModes + m + 1;
            nInd = numModes + n + 1;
            dx = 1i*(n + f);
            if abs(m-n)<=numModes
                mnInd = N + 1 + m - n;
                RR21(mInd,nInd) = coshTerm1Hat(mnInd);
                LL21(mInd,nInd) = -1i*k*coshTerm1Hat2(mnInd);
                LL22(mInd,nInd) = -1i*dx/mu*sinhTerm1Hat(mnInd);
                
                RR31(mInd,nInd) = coshTerm2Hat(mnInd);
                LL31(mInd,nInd) = -1i*k*coshTerm2Hat2(mnInd);
                LL33(mInd,nInd) = -1i*dx/mu*sinhTerm2Hat(mnInd);
            end
        end
        
    end
    
    
    LL11 = L11 + L11m.*DX;
    LL12 = L12m.*DX;
    LL13 = L13m.*DX;
    RR11 = R11om./DX + R11;
    RR12 = R12;
    RR13 = R13;
    
    if f==0
        selectModes = [-numModes:-1 1:numModes] + numModes + 1;
        RR11 = RR11(selectModes,selectModes);
        RR12 = RR12(selectModes,selectModes);
        RR13 = RR13(selectModes,selectModes);
        
        RR21 = RR21(selectModes,selectModes);
        RR31 = RR31(selectModes,selectModes);
        LL11 = LL11(selectModes,selectModes);
        LL12 = LL12(selectModes,selectModes);
        LL13 = LL13(selectModes,selectModes);
        LL21 = LL21(selectModes,selectModes);
        LL22 = LL22(selectModes,selectModes);
        LL31 = LL31(selectModes,selectModes);
        LL33 = LL33(selectModes,selectModes);
        
    end
    
    A0 = (-LL11 + LL12*inv(LL22)*LL21 + LL13*inv(LL33)*LL31);
    A1 = (RR11 - RR12*inv(LL22)*LL21 - LL12*inv(LL22)*RR21 - LL13*inv(LL33)*RR31 - RR13*inv(LL33)*LL31);
    A2 =  (RR12*inv(LL22)*RR21 + RR13*inv(LL33)*RR31);
    
    [eVect,eVal,SS] = polyeig(A0,A1,A2);
    if f==0
        eVals(:,fInd) = [0;0;eVal];
    else
        eVals(:,fInd) = (eVal);
        
    end
    
end