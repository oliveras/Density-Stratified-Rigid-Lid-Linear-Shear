function integral = integrate(f,dx)
%
% 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% I = integrate(f,dx)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Assuming that f is a periodic function with one but not both endpoints 
% not included in the range, the function 'integrate' returns the integral
% calculated over one period by appending the first value of the function
% to the end of the array and then using the trapezoid rule.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Katie Oliveras 4/23/2016


g = [f;f(1)];
integral = dx/2*(2*sum(g) - g(1) - g(end));

integral = dx/2*(2*sum(f));