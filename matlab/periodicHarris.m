
% periodic hyperbolic tangent profile

% hyperbolic tangent periodic profile
%
Ly=4*pi;
w0=.5;
w0inv=1/w0;
y=[0:.01:1]*Ly;
f = tanh((y-.25*Ly)*w0inv) - ...
    tanh((y-.75*Ly)*w0inv) - 1;
F = w0*log(cosh((y-.25*Ly)*w0inv)) - ...
    w0*log(cosh((y-.75*Ly)*w0inv)) - y;
% this loop does essentially nothing; theoretically,
% k=1:infty makes the function smooth at the boundaries.
for k=1:0
  for j=[-k,k]
    f=f+tanh((y-(j+.5)*Ly).*w0inv)-tanh((y-(j+.75)*Ly).*w0inv);
  F=F+w0*log(cosh((y-(j+.25)*Ly)*w0inv)) - ...
      w0*log(cosh((y-(j+.75)*Ly)*w0inv));
  end
end


figure(2);
plot(y,f);
figure(3);
plot(y,F);

% periodic version of the hyperbolic tangent
% for use in defining a doubly periodic harris sheet
%function out=tanhp(p,x)
%  out = tanh(x);
%  n = 4; % n approximates infinity and should be chosen based on p
%  % This series is only valid in the interval [-p,p]
%  for j=1:n
%    % need to compute this in a different way
%    % to avoid loss of precision
%    correction = tanh(x-2*j*p) - tanh(x-(2*j+1)*p) ...
%    out = out - correction;
%  end
%  %correction = tanh(x+p) + tanh(x-p);
%  %out = out - correction;
%end

