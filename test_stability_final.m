% Computes the stability of a small amplitude traveling wave solution.

clear all
close all
clc


addpath ./subroutines -end

% ===========================================================================================================
% Model parameters
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% set variable values here
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
N = 8;                     M = 2*N + 1;
L = 2*pi;                   x = (0:2*N)'/M*L;
delta = 2;                  rho = 1.028;
epsilon = 1;                w1 = -1;
mu = 2;                     w2 = 1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% ===========================================================================================================
% Initializing the guess for the solution (use continuation methods to
% calculate all but weakly nonlinear solutions).
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Initializing the bifurcation branch
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
k0 = 1;  amp = .05; 

T = tanh(mu*k0)*tanh(mu*k0*delta)/(mu*k0*(rho*tanh(mu*k0)+tanh(mu*k0*delta)));
c0 = 1/2*(-T*(w1-rho*w2) + sqrt((w1-rho*w2)^2*T^2 + 4*(rho-1)*T));
eta0 = amp*cos(k0*x); q2x0 = (mu*k0*c0/tanh(mu*k0*delta) + w2)*eta0; iniGuessX = [hat(eta0); hat(q2x0); c0];
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
p.amp = amp;        p.qSign=qSign;     
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ===========================================================================================================


% ===========================================================================================================
% Items important and relative to the numerical construction of solutions
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Solver parameters
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tol = 1e-15;
options = optimoptions('fsolve','Display','off','Jacobian','on','TolX',tol,'TolFun',tol,'Algorithm','levenberg-marquardt');

p.ContinuationType = 2;  % 1 = FourierMode, 2 = InfNorm, 3 = Max Value

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ===========================================================================================================



% ===========================================================================================================
% Compute the solution
% ===========================================================================================================
[X,fval,exitflag,output,jacobian] = fsolve(@(X) getSpectralIntegralsNonDimensional(X,p,options),iniGuessX,options);
% ===========================================================================================================


% ===========================================================================================================
% Evaluate the spectra
% ===========================================================================================================
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Choose the Floquet Parameters and the truncation of the eigenvalue system
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% NOTE: p.numModes must be less than or equal to the number of modes used
% for the solution.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
floquetParams = linspace(0,.5,500); p.numModes = N;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Calulate and Plot the Spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[eVals] = getSpectralMatricesNonDimensionalParFor(X,p,floquetParams);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot the Spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subplot(2,1,1)
plot(floquetParams,abs(real(eVals)),'.k');
xlabel('Floquet Parameters - $p$','Interpreter','Latex','Fontsize',14)
ylabel('Growth Rate - Re$(\lambda)$','Interpreter','Latex','Fontsize',14)
set(gca,'TickLabelInterpreter','Latex','FontSize',12)

subplot(2,1,2)
plot(eVals,'.k')
hold on
plot(conj(eVals),'.k')
xlabel('Re$(\lambda)$','Interpreter','Latex','Fontsize',14)
ylabel('Im$(\lambda)$','Interpreter','Latex','Fontsize',14)
set(gca,'TickLabelInterpreter','Latex','FontSize',12)

axis tight
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% ===========================================================================================================
