function [etaOut, q1xOut, q2xOut,cOut, cRatio] = plotNonDimensionalSolutions(X,p,ampIter,plotFig,iniGuessX,muVect,eVals)
plotsOn=1;
% ===========================================================================================================
% Formatting for output to screen
% ===========================================================================================================
linespace = '=========================================================================================================';
Linespace = '---------------------------------------------------------------------------------------------------------';
% ===========================================================================================================

% ===========================================================================================================
% Pulling the parameter values from the structured array
% ===========================================================================================================
N = p.N; M = p.M; amp = p.amp; w1 = p.w1; w2 = p.w2; rho = p.rho;
ep = p.epsilon; mu = p.mu; delta = p.delta; qSign=p.qSign; numSols = p.numSols;
x = (0:2*N)'/M*p.L;
% ===========================================================================================================

% ===========================================================================================================
% Initial Guesses for plotting comparisons
% ===========================================================================================================
etaHat0 = iniGuessX(1:M); q2x0Hat = iniGuessX(M + (1:M)); c0 = iniGuessX(end);
eta0 = (invHat(etaHat0));
% ===========================================================================================================


% ===========================================================================================================
% Final solutions for plotting comparisons
% ===========================================================================================================
etaHat = X(1:M); q2xHat = X(M + (1:M)); c = X(end);
eta = real(invHat(etaHat)); etaX = d(eta);
q2x = real(invHat(q2xHat)); Q2 = ep*q2x - w2*ep*eta - c;
q1x=(w1*eta*ep + c + qSign*sqrt(rho*(Q2).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + (ep*mu*etaX).^2)));
etaHat = hat(eta); q2xHat = hat(q2x);
% ===========================================================================================================

% ===========================================================================================================
% Testing to make sure that solution was indeed real valued.  Warning printed if a wave-speed c was detected
% with imaginary part greater than 10^-10
if abs(imag(c))>1e-10
    disp('complex c')
else
    c = real(c);
end
% ===========================================================================================================

if plotsOn
    % ===========================================================================================================
    % ===========================================================================================================
    % Creating Plots
    % ===========================================================================================================
    %------------------------------------------------------------------------------------------------------------
    % Plotting the solution with two fluids visualized
    %------------------------------------------------------------------------------------------------------------
    subplot(2,2,1:2)
    %------------------------------------------------------------------------------------------------------------
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Creating the filled polygons for the fluid layers
    x2 = [x; 2*pi; 2*pi; flipud(x)];
    bottomY = [(min(eta)-delta/(delta+1)*amp)*ones(M+1,1); eta(1); flipud(eta)];
    topY = [(max(eta)+delta/(delta+1)*amp)*ones(M+1,1); eta(1); flipud(eta)];
    fill(x2,bottomY,[0 0 .8],'FaceAlpha',.4); hold on
    fill(x2,topY,[0 0 .8],'FaceAlpha',.1)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plot([x;2*pi],[eta;eta(1)],'k','LineWidth',2);
    plot([x;2*pi],[amp*cos(x);amp*cos(x(1))],'r--','LineWidth',1.5); hold off
    subplot(2,2,1:2)
    %     [Xx Zz psiT psiB pT pB pTe pBe x_hr eta_hr] = getPsi(X,p,[-p.delta,1],500);
    %     psiB = psiB - min(min(psiB)); psiB = psiB/max(max(psiB));
    %     psiT = psiT - min(min(psiT)); psiT = psiT/max(max(psiT));
    %     psiT = real(psiT); psiB = real(psiB);
    %     contour(Xx,Zz,psiT,1-cos(linspace(0,pi/2,10)),'k')
    %     hold on
    %     contour(Xx,Zz,psiB,sin(linspace(0,pi/2,10)),'k')
    %     plot(x_hr,eta_hr,'k','LineWidth',3)
    
    xlabel('$x$','FontSize',18,'Interpreter','Latex')
    ylabel('$z$','FontSize',18,'Interpreter','Latex')
    title(['Solution with $\omega_1 = $ ',num2str(w1),' and $\omega_2 = $ ',num2str(w2)],'Interpreter','Latex','FontSize',18)
    
    axis tight
    xlabel('$x$','FontSize',18,'Interpreter','Latex'); ylabel('$z$','FontSize',18,'Interpreter','Latex')
    title(['Solution with $\omega_1 = $ ',num2str(w1),' and $\omega_2 = $ ',num2str(w2)],'Interpreter','Latex','FontSize',18)
    
    hold off
    
    %------------------------------------------------------------------------------------------------------------
    
    %------------------------------------------------------------------------------------------------------------
    % Plotting the bifurcation curves of c vs different definitions of amplitude
    %------------------------------------------------------------------------------------------------------------
    subplot(2,2,3)
    %------------------------------------------------------------------------------------------------------------
    
    plot(c0,amp,'*r'); hold on
    plot(c,amp,'bo','MarkerFaceColor','b')
    if ampIter==1
        plot(c0,0,'bo');
    end
    
    %------------------------------------------------------------------------------------------------------------
    % Taking a close look at the Fourier Coefficients of the solution for eta to ensure convergence
    %------------------------------------------------------------------------------------------------------------
    subplot(2,2,4)
    %------------------------------------------------------------------------------------------------------------
    if nargin<6
        stem(1:N,real(etaHat(N+2:end))); hold on
        stem(1:N,real(etaHat0(N+2:end))); hold off
        axis tight
        ylim([-1 1]*1e-10)
    else
        plot(muVect,real(eVals),'xk');
        axis tight
    end
    pause(.1)
    %------------------------------------------------------------------------------------------------------------
    %------------------------------------------------------------------------------------------------------------
    
end
% ===========================================================================================================
% Saving Output
% ===========================================================================================================
Q2 = real(ep*q2x - w2*ep*eta - c); etaOut = real(eta); etaX = d(eta);
q1xOut = 1/ep*(w1*eta*ep + c + qSign*sqrt(rho*(Q2).^2 - (rho-1)*(c^2 - 2*ep*eta).*(1 + ep^2*mu^2*etaX.^2)));q1x = q1xOut;
Q1 = ep*q1x - w1*eta*ep - c; q2xOut = real(q2x);
cOut = c; maxEta = max(eta); minEta = min(eta);
% ===========================================================================================================

% ===========================================================================================================
% Calculating values that may be related to wave of maximum height
% ===========================================================================================================
if (p.topStagnation==1)&&(p.bottomStagnation==0)
    cRatio = max((1 - 1/rho)*(c^2 - 2*ep*eta).*(1 + (d(eta)*ep*mu).^2)./(ep*q2x - ep*w2*eta - c).^2);
elseif (p.topStagnation==0)&&(p.bottomStagnation==1)
    cRatio = max(ep*w2*eta./(ep*q2x - c));
else
    cRatio = min(max((1 - 1/rho)*(c^2 - 2*ep*eta).*(1 + (d(eta)*ep*mu).^2)./((ep*q2x - ep*w2*eta - c).^2)),max(ep*w2*eta./(ep*q2x - c)));
end
cRatio = ([max((1 - 1/rho)*(c^2 - 2*ep*eta).*(1 + (d(eta)*ep*mu).^2)./((ep*q2x - ep*w2*eta - c).^2)),max(ep*w2*eta./(ep*q2x - c))]);
% ===========================================================================================================


% ===========================================================================================================
% Displaying output to screen in terms of raw data
% ===========================================================================================================
if ampIter==1
    disp(linespace); disp(linespace);
    strOut = sprintf('  Iter # / %-5d    |   Amplitude   |   Infinity Norm    |         Ratio C         |      Residual     |',numSols);
    disp(strOut);
    disp(Linespace);
end
strOut = sprintf('  k = %5d / %-5d | Amp = %0.5f | Inf Norm = %0.5f | C Ratio = %-3.4f,%-3.4f | Error = %1.3e | ',ampIter,numSols,1/2*sum((ones(M,1) - (-1).^((-N:N)')).*etaHat),max(eta),cRatio(1), cRatio(2),norm(getSpectralIntegralsNonDimensional(X,p)));
cRatio = min(cRatio);
disp(strOut)
% ===========================================================================================================

