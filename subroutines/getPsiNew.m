function [Xx Zz psiT psiB pT pB x_hr eta_hr pressureEtaT pressureEtaB topmask botmask pbot_hr] = getPsiNew(X,p,zRange,res)
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



sumTermT = 0; sumTermTX = 0; sumTermTZ = 0;
sumTermB = 0; sumTermBX = 0; sumTermBZ = 0;
psizb = 0;
for k = -N:N
    if abs(k)>0
        sinhTermT = sinh(mu*k*(ep*eta - 1));    coshTermT = cosh(mu*k*(ep*eta - 1));
        sinhTermB = sinh(mu*k*(ep*eta + del));  coshTermB = cosh(mu*k*(ep*eta + del));
        
        psiHatT(k+N+1) = hat(1/(mu*k*cosh(mu*k))*(w1/(mu*k)*sinh(mu*k*ep*eta) + sqrtTerm.*cosh(mu*k*ep*eta)),k);
        psiHatB(k+N+1) = hat(1/(mu*k*cosh(mu*k*del))*(w2/(mu*k)*sinh(mu*k*ep*eta) + (ep*q2x-w2*eta-c).*cosh(mu*k*ep*eta)),k);
        
        sumTermT = sumTermT + psiHatT(k+N+1)*exp(1i*k*(x)).*sinh(mu*k*(ep*eta - 1));
        sumTermTX = sumTermTX + (1i*k)*psiHatT(k+N+1)*exp(1i*k*(x)).*sinh(mu*k*(ep*eta - 1));
        sumTermTZ = sumTermTZ + (mu*k)*psiHatT(k+N+1)*exp(1i*k*(x)).*cosh(mu*k*(ep*eta - 1));
        
        sumTermB = sumTermB + psiHatB(k+N+1)*exp(1i*k*(x)).*sinh(mu*k*(ep*eta + del));
        sumTermBX = sumTermBX + (1i*k)*psiHatB(k+N+1)*exp(1i*k*(x)).*sinh(mu*k*(ep*eta + del));
        sumTermBZ = sumTermBZ + (mu*k)*psiHatB(k+N+1)*exp(1i*k*(x)).*cosh(mu*k*(ep*eta + del));
        
        psizb = psizb + psiHatB(k+N+1)*mu*k*exp(1i*k*x);
        
    else
    a1T = hat(w1*(ep*eta-1) + sqrtTerm,0);
    a1B = hat(w2*(ep*eta+del) + (ep*q2x - w2*ep*eta - c),0);
    end
        
end


sumTermT = real(sumTermT);  sumTermTX = real(sumTermTX);    sumTermTZ = real(sumTermTZ);
sumTermB = real(sumTermB);  sumTermBX = real(sumTermBX);    sumTermBZ = real(sumTermBZ);

a2T = -w1/2; a2B = -w2/2;
a0T = - hat(a1T*(eta-1) + a2T*(eta-1).^2 + sumTermT,0);

a0B = - hat(a1B*(eta+del) + a2B*(eta+del).^2 + sumTermB,0);


psizb = a1B+real(psizb);

pbot = 1/2*(c^2 + 2*del - 2*w2*a0B - psizb.^2);
psiEtaT = a0T + a1T*(eta-1) + a2T*(eta-1).^2 + sumTermT;
psiEtaB = a0B + a1B*(eta+del) + a2B*(eta+del).^2 + sumTermB;

psiTX = sumTermTX;
psiTZ = a1T + a2T*2*(eta-1)+sumTermTZ;

psiBX = sumTermBX;
psiBZ = a1B + a2B*2*(eta+del)+sumTermBZ;

pressureEtaT = 1/2*c^2 - 1/2*(mu^2*psiTX.^2 + psiTZ.^2) - ep*eta;
pressureEtaB = rho*(1/2*c^2 - 1/2*(mu^2*psiBX.^2 + psiBZ.^2) - ep*eta);



%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------------
if nargin<3
    x_hr = linspace(0,4*pi,501); z_hr = linspace(-del,1,200);
    
else
    x_hr = linspace(0,4*pi,res);
    z_hr = linspace(zRange(1), zRange(2),res);
end
eta_hr = 0; pbot_hr = 0; pHat = hat(pbot);
for n = -N:N
    pbot_hr = pbot_hr + exp(1i*n*x_hr)*(pHat(n+N+1));
    eta_hr = eta_hr + exp(1i*n*x_hr)*(etaHat(n+N+1));
end
eta_hr = real(eta_hr);
pbot_hr = real(pbot_hr);


[Xx Zz] = meshgrid(x_hr,z_hr);
[topmask, botmask] = createMasksNonDimensional(Zz,eta_hr,ep,del);


sumTermT = 0; sumTermTX = 0; sumTermTZ = 0;
sumTermB = 0; sumTermBX = 0; sumTermBZ = 0;
for k = -N:N
    if abs(k)>0
        sumTermT = sumTermT + psiHatT(k+N+1)*exp(1i*k*(Xx)).*sinh(mu*k*(Zz - 1));
        sumTermTX = sumTermTX + (1i*k)*psiHatT(k+N+1)*exp(1i*k*(Xx)).*sinh(mu*k*(Zz - 1));
        sumTermTZ = sumTermTZ + (mu*k)*psiHatT(k+N+1)*exp(1i*k*(Xx)).*cosh(mu*k*(Zz - 1));
        
        sumTermB = sumTermB + psiHatB(k+N+1)*exp(1i*k*(Xx)).*sinh(mu*k*(Zz + del));
        sumTermBX = sumTermBX + (1i*k)*psiHatB(k+N+1)*exp(1i*k*(Xx)).*sinh(mu*k*(Zz + del));
        sumTermBZ = sumTermBZ + (mu*k)*psiHatB(k+N+1)*exp(1i*k*(Xx)).*cosh(mu*k*(Zz + del));
  
    end
        
end

sumTermT = real(sumTermT);  sumTermTX = real(sumTermTX);    sumTermTZ = real(sumTermTZ);
sumTermB = real(sumTermB);  sumTermBX = real(sumTermBX);    sumTermBZ = real(sumTermBZ);


psiT = a0T + a1T*(Zz-1) + a2T*(Zz-1).^2 + sumTermT;
psiB = a0B + a1B*(Zz+del) + a2B*(Zz+del).^2 + sumTermB;

psiTX = sumTermTX;
psiTZ = a1T + a2T*2*(Zz-1)+sumTermTZ;

psiBX = sumTermBX;
psiBZ = a1B + a2B*2*(Zz+del)+sumTermBZ;

pT = 1/2*c^2 - 1/2*(mu^2*psiTX.^2 + psiTZ.^2) - Zz - w1*psiT;
pB = rho*(1/2*c^2 - 1/2*(mu^2*psiBX.^2 + psiBZ.^2) - Zz - w2*psiB);


%psiT = psiT.*topmask; pT = pT.*topmask;
%psiB = psiB.*botmask; pB = pB.*botmask;
