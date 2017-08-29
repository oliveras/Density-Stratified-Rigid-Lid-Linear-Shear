function [jacobian] = getJacobianNonDimensional(X,p)
N = p.N; M = p.M; amp = p.amp;
x = (0:2*N)'/M*2*pi; dx = x(2)-x(1);
w1 = p.w1; w2 = p.w2; rho = p.rho;
ep = p.epsilon; mu = p.mu; del = p.delta;
qSign=p.qSign;

Q1XSign=p.qSign;
%
etaHat = X(1:M); q2xHat = X(M+(1:M)); c = (X(end));
%
eta = invHat(etaHat);
etaX = d(eta);
q2x = invHat(q2xHat);
Q2 = ep*q2x - w2*ep*eta - c - p.U0; 
dQ2dq2 = ep;
Q1 = qSign* sqrt(rho*(Q2).^2 + (1-rho)*(c^2 - 2*ep*eta).*(1 + (ep*mu)^2*etaX.^2));
etaX = d(eta);

dQ2dc = -1;
dQ1dc = 1./Q1.*(rho*Q2.*dQ2dc - (rho-1)*c*(1 + ep^2*mu^2*etaX.^2));
dQ2dEta = -w2*ep;
dQ1dq2 = 1./Q1.*(rho*Q2)*ep;

jacobian = zeros(2*M+1,2*M+1);

for n = -N:N
    nInd = n + N + 1;
    k = n;
    coshTerm1 = cosh(mu*k*ep*eta) - sinh(mu*k*ep*eta)*tanh(mu*k);
    sinhTerm1 = sinh(mu*k*ep*eta) - cosh(mu*k*ep*eta)*tanh(mu*k);
    coshTerm2 = cosh(mu*k*ep*eta) + tanh(mu*k*del)*sinh(mu*k*ep*eta);
    sinhTerm2 = sinh(mu*k*ep*eta) + tanh(mu*k*del)*cosh(mu*k*ep*eta);
    
    
    if (abs(n)>0)
        for m = -N:N
            mInd = N + m + 1;
            l = m;
            if abs(n-m)<=N
                dQ1dEta = 1./Q1.*(rho*Q2*dQ2dEta - (rho-1)*(1i*l*etaX.*(c^2-2*ep*eta)*ep^2*mu^2 - ep*(1 + ep^2*mu^2*etaX.^2)));
                
                dAdEta = w1*sinhTerm1*ep + (mu*k)*coshTerm1.*Q1*ep + sinhTerm1.*dQ1dEta;
                dBdEta = w2*sinhTerm2*ep + (mu*k)*coshTerm2.*Q2*ep + sinhTerm2.*dQ2dEta;
                dAdq = sinhTerm1.*dQ1dq2;
                dBdq = sinhTerm2.*dQ2dq2;
                jacobian(nInd,mInd) = integrate(exp(1i*(l-k)*x).*dAdEta,dx);%hat(dAdEta,n-m);%
                jacobian(nInd+M,mInd) = integrate(exp(1i*(l-k)*x).*dBdEta,dx);%hat(dBdEta,n-m);%
                jacobian(nInd,mInd+M) = integrate(exp(1i*(l-k)*x).*dAdq,dx);%hat(dAdq,n-m);%
                jacobian(nInd+M,mInd+M) = integrate(exp(1i*(l-k)*x).*dBdq,dx);%hat(dBdq,n-m);%
%                 if n>0
%                     jacobian(nInd,mInd) = real(jacobian(nInd,mInd));
%                     jacobian(nInd+M,mInd) = real(jacobian(nInd+M,mInd));
%                     jacobian(nInd,mInd+M) = real(jacobian(nInd,mInd+M));
%                     jacobian(nInd+M,mInd+M) = real(jacobian(nInd+M,mInd+M));
%                 else
%                     jacobian(nInd,mInd) = imag(jacobian(nInd,mInd));
%                     jacobian(nInd+M,mInd) = imag(jacobian(nInd+M,mInd));
%                     jacobian(nInd,mInd+M) = imag(jacobian(nInd,mInd+M));
%                     jacobian(nInd+M,mInd+M) = imag(jacobian(nInd+M,mInd+M));
%                 end
                
            end
        end
    end
    jacobian(nInd,end) = integrate(exp(-1i*k*x).*sinhTerm1.*dQ1dc,dx);%hat(sinhTerm1.*dQ1dc,n);%
    jacobian(nInd+M,end) = integrate(exp(-1i*k*x).*sinhTerm2.*dQ2dc,dx);%hat(sinhTerm2.*dQ2dc,n);%
%     if n>0
%         jacobian(nInd,end) = real(jacobian(nInd,end));
%         jacobian(nInd+M,end) = real(jacobian(nInd,end));
%     else
%         jacobian(nInd,end) = imag(jacobian(nInd,end));
%         jacobian(nInd+M,end) = imag(jacobian(nInd+N,end));
%     end
    
end



jacobian(N+1,N+1) = 1;
jacobian(M+N+1,M+N+1) = 1;%[-w2*ones(1,M) ones(1,M) -1];
jacobian(end,1:M) = .5*(ones(1,M) - (-1).^((-N:N)));


