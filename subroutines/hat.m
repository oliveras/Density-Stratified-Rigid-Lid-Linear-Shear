function uHat = hat(u,k);
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% uHat = hat(u,k)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% Determines the Fourier coefficients of the Fourier series for the 
% periodic function 'u'.  The argument 'k' is optional.  When 'k' is
% defined, uHat only returns the Fourier coefficient corresponding to the
% 'exp(i*k*x)' term.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 4/23/2016


M = length(u);
N = (M-1)/2;
x = (0:M-1)'/M*2*pi;
dx = x(2);

uHat = ((fftshift(fft(u))/length(u)));

if nargin ==2
    k+N+1;
    uHat = uHat(k + N + 1);
end



