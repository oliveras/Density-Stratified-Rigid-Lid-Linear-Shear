function [Xx Zz psiT psiB pT pB pTe pBe x eta psiTUM psiBUM] = getPsi(X,p,zRange,res)
N = p.N; M = p.M;
c = X(end);
c = real(c);
eta = invHat(X(1:M));etaX = d(eta);
q2x = invHat(X(M+1:end-1));
x = (0:2*N)'/M*2*pi;
etaHat = hat(eta); q2xHat = hat(q2x);

rho = p.rho; ep = p.epsilon; del = p.delta;
mu = p.mu; w1 = p.w1; w2 = p.w2;
sqrtTerm = p.qSign*sqrt(rho*(ep*q2x - ep*w2*eta - c).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + ep^2*mu^2*etaX.^2));

kv = (-N:N)'; kv(N+1) = 1;


for k = -N:N
    if abs(k)>0
        sinhTermT = sinh(mu*k*(ep*eta - 1));    coshTermT = cosh(mu*k*(ep*eta - 1));
        sinhTermB = sinh(mu*k*(ep*eta)) + cosh(mu*k*ep*eta)*tanh(mu*k*del);
        coshTermB = cosh(mu*k*ep*eta) + sinh(mu*k*ep*eta)*tanh(mu*k*del);
        %sinhTermB = sinh(mu*k*(ep*eta + del)); coshTermB = cosh(mu*k*(ep*eta + del));
        
        %fHatT(N+k+1,1) = hat(w1*sinhTermT/(k*mu) + coshTermT.*sqrtTerm,k)/(mu*k*ep);
        %fHatB(N+k+1,1) = hat(w2*sinhTermB/(k*mu) + coshTermB.*(ep*q2x - ep*w2*eta - c),k)/(mu*k*ep);
        fHatT(N+k+1,1) = -1/(mu*k*sinh(mu*k))*hat(-sinh(mu*k*ep*eta).*(sqrtTerm)+w1/(mu*k)*cosh(mu*k*ep*eta),k);
        fHatB(N+k+1,1) = -1/(mu*k*cosh(mu*k*(-del)))*hat(-cosh(mu*k*ep*eta).*(ep*q2x - ep*w2*eta - c)-w2/(mu*k)*sinh(mu*k*ep*eta),k);
    end
end


eta = 0; q2x = 0; etaX = 0;

x = linspace(0,4*pi,res);
z = linspace(zRange(1),zRange(2),res);
for k = -N:N
    eta = eta + exp(1i*k*x)*etaHat(k + N + 1);
    etaX = etaX + 1i*k*exp(1i*k*x)*etaHat(k + N + 1);
    q2x = q2x + exp(1i*k*x)*q2xHat(k + N + 1);
end
eta = real(eta); etaX = real(etaX); q2x = real(q2x);

sqrtTerm = p.qSign*sqrt(rho*(ep*q2x - ep*w2*eta - c).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + ep^2*mu^2*etaX.^2));


[Xx Zz] = meshgrid(x,z);[topmask, botmask] = createMasksNonDimensional(Zz,eta,ep,del);

psiT = -c*Zz - 1/2*w1*Zz.^2; psiTx = 0; psiTz = -c - w1*Zz;
psiB = -c*Zz - 1/2*w2*Zz.^2; psiBx = 0; psiBz = -c - w2*Zz;
psiTe = -c*ep*eta - 1/2*w1*ep^2*eta.^2; psiTxe = 0; psiTze = -c - w1*ep*eta;
psiBe = -c*ep*eta - 1/2*w2*ep^2*eta.^2; psiBxe = 0; psiBze = -c - w2*ep*eta;

for k = -N:N
    if abs(k)>0
        psiT = psiT + fHatT(N+k+1)*exp(1i*k*Xx).*sinh(mu*k*(Zz - 1));
        psiTx = psiTx + fHatT(N+k+1)*1i*k*exp(1i*k*Xx).*sinh(mu*k*(Zz - 1));
        psiTz = psiTz + fHatT(N+k+1)*mu*k*exp(1i*k*Xx).*cosh(mu*k*(Zz - 1));
        
        psiB = psiB + fHatB(N+k+1)*exp(1i*k*Xx).*sinh(mu*k*(Zz + del));
        psiBx = psiBx + fHatB(N+k+1)*1i*k*exp(1i*k*Xx).*sinh(mu*k*(Zz + del));
        psiBz = psiBz + fHatB(N+k+1)*mu*k*exp(1i*k*Xx).*cosh(mu*k*(Zz + del));
        
        psiTe = psiTe + fHatT(N+k+1)*exp(1i*k*x).*sinh(mu*k*(ep*eta - 1));
        psiTxe = psiTxe + fHatT(N+k+1)*1i*k*exp(1i*k*x).*sinh(mu*k*(ep*eta - 1));
        psiTze = psiTze + fHatT(N+k+1)*mu*k*exp(1i*k*x).*cosh(mu*k*(ep*eta - 1));
        
        psiBe = psiBe + fHatB(N+k+1)*exp(1i*k*x).*sinh(mu*k*(ep*eta + del));
        psiBxe = psiBxe + fHatB(N+k+1)*1i*k*exp(1i*k*x).*sinh(mu*k*(ep*eta + del));
        psiBze = psiBze + fHatB(N+k+1)*mu*k*exp(1i*k*x).*cosh(mu*k*(ep*eta + del));
        
    end
end

if max(max(abs(imag(psiT.*topmask))))<1e-10
    psiT = real(psiT);
else
    disp('imaginary psi_top')
end
if max(max(abs(imag(psiB))))<1e-10
    psiB = real(psiB);
else
    disp('imaginary psi_bottom')
end
psiBUM = psiB; psiTUM = psiT;


pT = real(-1/2*(psiTx.^2*mu^2 + psiTz.^2) - w1*psiT - Zz);
pB = rho*real(-1/2*(psiBx.^2*mu^2 + psiBz.^2) - w2*psiB - Zz);


psiT = psiT.*topmask; pT = pT.*topmask;


pTe = real(-1/2*(psiTxe.^2*mu^2 + psiTze.^2) - w1*psiTe - ep*eta);
pTe = pTe - 1/(4*pi)*trapz(x,pTe);
pT = pT - 1/(4*pi)*trapz(x,pTe);

pBe = rho*real(-1/2*(psiBxe.^2*mu^2 + psiBze.^2) - w2*psiBe - ep*eta);
pBe = pBe - 1/(4*pi)*trapz(x,pBe);
pB = pB - 1/(4*pi)*trapz(x,pBe);

psiB = psiB.*botmask; pB = pB.*botmask;
