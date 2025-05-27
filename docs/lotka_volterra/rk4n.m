function [tout, yout] = rk4n(FunFcn, tspan, y0, N)

% RK4	Integrates a system of ordinary differential equations using
%	the fourth order Runge-Kutta  method.  See also ODE45 and
%	ODEDEMO.M.
%	[t,y] = rk4('yprime', tspan, y0) integrates the system
%	of ordinary differential equations described by the M-file
%	yprime.m over the interval tspan=[t0,tfinal] and using initial
%	conditions y0.
%	[t, y] = rk4(F, tspan, y0, N) uses number of points
%   Modified by J.bastien, Mai 2025
%
% INPUT:
% F     - String containing name of user-supplied problem description.
%         Call: yprime = fun(t,y) where F = 'fun'.
%         t      - Time (scalar).
%         y      - Solution column-vector.
%         yprime - Returned derivative column-vector; yprime(i) = dy(i)/dt.
% tspan = [t0, tfinal], where t0 is the initial value of t, and tfinal is
%         the final value of t.
% y0    - Initial value column-vector.
% ssize - The number of points. (Default: N = 100).
%
% OUTPUT:
% t  - Returned integration time points (column-vector).
% y  - Returned solution, one solution column-vector per tout-value.
%
% The result can be displayed by: plot(t,y).


% Initialization

t0=tspan(1);
tfinal=tspan(2);
%pm = sign(tfinal - t0);  % Which way are we computing?
if nargin < 4, N = 100; end
% Inutile
% if ssize < 0, ssize = -ssize; end
% h = pm*ssize;
h=(tfinal-t0)/(N-1);
%t = t0;
y = y0(:);

% We need to compute the number of steps.
% Inutile 
% dt = abs(tfinal - t0);
% N = floor(dt/ssize) + 1;
% if (N-1)*ssize < dt
%   N = N + 1;
% end

% Initialize the output.

tout = (linspace(t0,tfinal,N)).';
%tout(1) = t;
yout = zeros(N,size(y,1));
yout(1,:) = y.';

% The main loop
for k=2:N
    % while (k < N)
    %   if pm*(t + h - tfinal) > 0
    %     h = tfinal - t;
    %     tout(k+1) = tfinal;
    %   else
    %     tout(k+1) = t0 +k*h;
    %   end
    %   k = k + 1;
    t = tout(k-1);

    % Compute the slopes
    s1 = feval(FunFcn, t, y); s1 = s1(:);
    s2 = feval(FunFcn, t + h/2, y + h*s1/2); s2=s2(:);
    s3 = feval(FunFcn, t + h/2, y + h*s2/2); s3=s3(:);
    s4 = feval(FunFcn, t + h, y + h*s3); s4=s4(:);
    y = y + h*(s1 + 2*s2 + 2*s3 +s4)/6;
    yout(k,:) = y.';
end