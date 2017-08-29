function u = invHat(uHat);
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% u = invHat(uHat)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% 
% Determines the function 'u' based on the Fourier coefficients given by
% the vector uHat.  Assumes that the vector has length 2N+1 where uHat(1)
% represents the -N-th Fourier coefficient, uHat(N+1) represents the 0th
% Fourier coefficient, and uHat(end) represents the N-th Fourier
% coefficient.
%
% Uses complex Fourier series.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Katie Oliveras 4/23/2016



% 
% M = length(uHat);
% N = (M-1)/2;
% x = (0:2*N)'/M*2*pi;
% u = zeros(size(uHat));
% 
% for j = 1:M
%     n = j - N - 1;
%     u = u + exp(i*n*x)*uHat(j);
% end


u = (ifft(ifftshift(uHat))*length(uHat));
if max(abs(imag(u)))<1e-10
    u = real(u);
end
