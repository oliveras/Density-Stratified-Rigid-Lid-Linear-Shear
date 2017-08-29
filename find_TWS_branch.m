clear all
close all
clc
format long
format compact

addpath ./subroutines -end
warning('off','MATLAB:datetime:NonstandardSystemTimeZone')

% ===========================================================================================================
% Model parameters
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% set variable values here
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = 8;                     M = 2*N + 1;
L = 2*pi;                   x = (0:2*N)'/M*L;
delta =1;                   rho = 1.028;
epsilon = 1;                w1 = 0;
mu = sqrt(.1);              w2 = 0;
               
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Amplitude Values
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ampVals = linspace(.0001,.001,50);
amp = ampVals(1);           numSols = length(ampVals);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create file name for the solutions to be saved as in terms of the model parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filePart = ['sol_epsilon_',regexprep(num2str(epsilon),'\.','_')];
filePart = [filePart,'_mu_',regexprep(num2str(mu),'\.','_')];
filePart = [filePart,'_delta_',regexprep(num2str(delta),'\.','_')];
filePart = [filePart,'_rho_',regexprep(regexprep(num2str(rho),'\.','_'),'\-','m_')];
filePart = [filePart,'_w1_',regexprep(regexprep(num2str(w1),'\.','_'),'\-','m_')];
filePart = [filePart,'_w2_',regexprep(regexprep(num2str(w2),'\.','_'),'\-','m_')];
fileStart = filePart;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ===========================================================================================================


% ===========================================================================================================
% Initializing Solution Branch
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Initializing the bifurcation branch
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
T = tanh(mu)*tanh(mu*delta)/(mu*(rho*tanh(mu)+tanh(mu*delta)));
c0 = 1/2*(-T*(w1-rho*w2) + sqrt((w1-rho*w2)^2*T^2 + 4*(rho-1)*T));
eta0 = amp*cos(x); q2x0 = (mu*c0/tanh(mu*delta) + w2)*eta0; iniGuessX = [hat(eta0); hat(q2x0); c0];
qSign = -sign(c0);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Determining if there is a stagnation point within each fluid as a function of model parameters for the
% trivial flow
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
p.topStagnation = (c0/-w1>0)&&(c0/-w1<1);
p.bottomStagnation = (c0/-w2 > - delta)&&(c0/-w2<0);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Saving model parameters into a structured array
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
p.U0 = 0;           p.N = N;            p.M = M;                p.x = x;
p.L = L;            p.c0 = c0;          p.epsilon = epsilon;    p.mu = mu;
p.delta = delta;    p.rho = rho;        p.w1 = w1;              p.w2 = w2;
p.amp = amp;        p.qSign=qSign;      p.numSols = numSols;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ===========================================================================================================


% ===========================================================================================================
% Items important and relative to the numerical construction of solutions
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Solver parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
NMAX = 64;          dN = 4;             tol = 1e-16;            histInd = 4;
options = optimoptions('fsolve','Display','off','Jacobian','on','TolX',tol,'TolFun',tol,'Algorithm','levenberg-marquardt');

startVal = 1;       XOut = zeros(2*(2*64+1)+1,numSols);  plotOn = 1; %Set to 0 to turn plots off
p.ContinuationType = 2;  %1 = FourierMode, 2 = InfNorm, 3 = Max Value

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ===========================================================================================================
% ===========================================================================================================


for ampIter = startVal:numSols
    
    if ampIter==1
        ampStart = ampVals(ampIter);
        plotFig = figure();
        iniGuessX(end) = iniGuessX(end);
    else
        amp = ampVals(ampIter);
        p.amp = amp;
        
        % Crude continuation method using splines to speed computation.
        if ampIter > histInd
            cGuess = spline(ampOut(ampIter-histInd:ampIter-1),cOut(ampIter-histInd:ampIter-1),amp);
      
        else
            if ampIter==2
                slope = ampOut(1)/(cOut(1)-p.c0);
                cGuess = (ampVals(2)-ampVals(1))/slope+cOut(1);
            else
                
                cGuess = spline([0 ampOut(1:ampIter-1)],[p.c0 cOut(1:ampIter-1)],amp);
                oldPoint = [cOut(ampIter-1) ampOut(ampIter-1)];
                newPoint = [cGuess amp]; dist = newPoint-oldPoint;
                slope = dist(2)/dist(1);
                distNorm = norm(dist);
                ampIncrement = ampVals(ampIter)-ampVals(ampIter-1);
                cGuess = cOut(ampIter-1)+ampIncrement*cos(atan2(dist(2),dist(1)));
            end
            
        end
        eta0 = invHat(X(1:M)); q2x0 = invHat(X(M+(1:M)));
        
        eta0 = amp*eta0/halfMaxMinOut(ampIter-1); 
        q20x = q2x0*amp/halfMaxMinOut(ampIter-1);
        eta0Hat = hat(eta0); q2x0Hat = hat(q2x0);
        
        
        if (norm(eta0Hat(M-dN/2:M))>1e-12)
            if N < NMAX
                Nold = N; N = Nold + dN; M = 2*N+1; x = (0:2*N)'/M*L;
                eta0Hat = [zeros(dN,1);eta0Hat;zeros(dN,1)];
                q2x0Hat = [zeros(dN,1);q2x0Hat;zeros(dN,1)];
                p.N = N; p.M = M; p.x = x;
            end
        end
        
        iniGuessX = [eta0Hat;q2x0Hat; cGuess];
        
    end
    iniGuessX = real(iniGuessX);
    
    [X,fval,exitflag,output,jacobian] = fsolve(@(X) getSpectralIntegralsNonDimensional(X,p,options),iniGuessX,options);
    if plotOn==1
        [eta0, q1x0, q2x0,c0]=plotNonDimensionalSolutions(X,p,ampIter,plotFig,iniGuessX);
        pause(.01)
    end
    
    XOut(1:length(X),ampIter) = X;
    NOut(ampIter) = N;
    cOut(ampIter) = real(X(end));
    eta = real(invHat(X(1:M)));
    halfMaxMinOut(ampIter) = .5*(max(eta)-min(eta));
    ampOut(ampIter) = real(1/2*sum((ones(M,1) - (-1).^((-N:N)')).*hat(eta)));
    

end
ampEnd = amp;
filePart = ['epsilon_',regexprep(num2str(epsilon),'\.','_')];
filePart = [filePart,'_mu_',regexprep(num2str(mu),'\.','_')];
filePart = [filePart,'_delta_',regexprep(num2str(delta),'\.','_')];
filePart = [filePart,'_w1_',regexprep(regexprep(num2str(w1),'\.','_'),'\-','m_')];
filePart = [filePart,'_w2_',regexprep(regexprep(num2str(w2),'\.','_'),'\-','m_')];
filePart = [filePart,'_rho_',regexprep(regexprep(num2str(rho),'\.','_'),'\-','m_')];
ampPart = ['_ampStart_',regexprep(regexprep(num2str(ampStart),'\.','_'),'\-','m_')];
ampPart = [ampPart,'_ampEnd_',regexprep(regexprep(num2str(ampEnd),'\.','_'),'\-','m_')];
fileName = ['./solutions/solutions_for_paper_',filePart,ampPart];
%return
disp(fileName)
save(fileName,'XOut','p','NOut','cOut','ampOut','halfMaxMinOut','filePart')
% selectedSols = round(numSols*[.25 .5 .85]); % Solutions to plot at the end
%plotNonDimensionalSolutionsForPaper(fileName,selectedSols);
