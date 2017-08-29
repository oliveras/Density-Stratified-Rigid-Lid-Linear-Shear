function du = d(u,L)

%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% u = d(u,L)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% Computes the derivative of the function u that is L-periodic.  If L is
% not specified, defaults to 2*pi periodic functions.  Uses the Fourier
% transforms and assumes that the function u is real-valued.  
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 4/23/2016

if nargin==1
    L = 2*pi;
end

uHat= hat(u); M = length(u); N = (M-1)/2;  kV = (-N:N)'*2*pi/L;

du = invHat(1i*kV.*uHat);
du = real(du);
