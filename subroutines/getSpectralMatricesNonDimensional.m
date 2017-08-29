function [eVals eVecs] = getSpectralMatricesNonDimensional(X,p,muVect)






N = p.N; M = p.M; amp = p.amp;
numModes = p.numModes;

w1 = p.w1; w2 = p.w2; rho = p.rho;
ep = p.epsilon; mu = p.mu; del = p.delta;
del1 = p.delta1;
qSign=p.qSign;
etaHat = X(1:M); q2xHat = X((M+1):(2*M)); c = (X(end));


eta = invHat(etaHat);
etaX = d(eta);
q2x = invHat(q2xHat);
Q2x = ep*q2x - ep*w2*eta-c;

Q1x = qSign*sqrt(rho*(Q2x).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + (ep*mu*etaX).^2));

etaX = d(eta); etaXX = d(etaX);
[L11 L11m L12m L13m R11om R11 R12 R13] = getConstantMatricesNonDimensional(X, p);

for muInd = 1:length(muVect)
    f = muVect(muInd); DX = 1i*(f + meshgrid(-numModes:numModes));
    
    for m = -numModes:numModes
        k = f + m;
        if p.depth==inf
            coshTerm1 = exp(-abs(k)*eta);
            sinhTerm1 = -sign(k)*exp(-abs(k)*eta);
            coshTerm2 = exp(abs(k)*eta);
            sinhTerm2 = sign(k)*exp(abs(k)*eta);
            
        else
            coshTerm1 = (cosh(mu*ep*k*eta) - tanh(mu*k*del1)*sinh(mu*ep*k*eta));
            sinhTerm1 = (sinh(mu*ep*k*eta) - tanh(k*mu*del1)*cosh(mu*ep*k*eta));
            coshTerm2 = (cosh(mu*ep*k*eta) + tanh(k*mu*del)*sinh(mu*ep*k*eta));
            sinhTerm2 = (sinh(mu*ep*k*eta) + tanh(k*mu*del)*cosh(mu*ep*k*eta));
            %coshTerm1 = cosh(mu*k*(ep*eta-1)); coshTerm2 = cosh(mu*k*(ep*eta+del));
            %sinhTerm1 = sinh(mu*k*(ep*eta-1)); sinhTerm2 = sinh(mu*k*(ep*eta+del));
            
        end
        
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
    ZM = zeros(size(RR11));
    R = [RR11 RR12 RR13; RR21 ZM ZM; RR31 ZM ZM];
    L = [LL11 LL12 LL13; LL21 LL22 ZM; LL31 ZM LL33];
    
    
    
    if f==0
        A0 = (-LL11 + LL12*inv(LL22)*LL21 + LL13*inv(LL33)*LL31);
        A1 = (RR11 - RR12*inv(LL22)*LL21 - LL12*inv(LL22)*RR21 - LL13*inv(LL33)*RR31 - RR13*inv(LL33)*LL31);
        A2 =  (RR12*inv(LL22)*RR21 + RR13*inv(LL33)*RR31);
        [eVect,eVal,SS] = polyeig(A0,A1,A2);
        
        eVals(:,muInd) = [0;0;eVal];
    else
        A0 = (-LL11 + LL12*inv(LL22)*LL21 + LL13*inv(LL33)*LL31);
        A1 = (RR11 - RR12*inv(LL22)*LL21 - LL12*pinv(LL22)*RR21 - LL13*inv(LL33)*RR31 - RR13*inv(LL33)*LL31);
        A2 =  (RR12*inv(LL22)*RR21 + RR13*inv(LL33)*RR31);
        
        %A0 = (-LL11 + LL12*(LL22)\LL21 + LL13*(LL33)\LL31);
        %A1 = (RR11 - RR12*(LL22)\LL21 - LL12*(LL22)\RR21 - LL13*(LL33)\RR31 - RR13*(LL33)\LL31);
        %A2 =  (RR12*(LL22)\RR21 + RR13*(LL33)\RR31);
        L = [A2 0*A2; 0*A2 eye(size(A2)) ];
        R = [-A1 -A0; eye(size(A2)) zeros(size(A2))];
        [eVect,eVal,SS] = polyeig(A0,A1,A2);   eVals(:,muInd) = (eVal);
        %[eVect, eVal] = eig( R,L,'chol');    eVals(:,muInd) = diag(eVal);
        
        
        eVecs(:,:,muInd) = eVect;
        
        
    end
    
end