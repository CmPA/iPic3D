% % INPUT
% % n: number of entries in the rainbow colormap
% % dx: it controls the amount of pure blue and red used at the 
% %     beginning and the end of the colormap, respectively

% % OUTPUT
% % cmapRainbow: Rainbow color map values of Red, Green, Blue
% % g: scalar values between dx &amp; (6-dx)
% % f: n values between 0 and 1. g is evaluated at these values

function [cmapRainbow, g, f]=rainbow_cmap(n, dx)

f = linspace(0,1,n); % generates n points between 0 and 1.

cmapRainbow=[];
for i=1:n
    g(i) = (6 - 2*dx) * f(i) + dx ; %scale f between dx and (6-dx)
    R = max(0, (3 - abs(g(i) - 4) - abs(g(i) - 5))/2); 
    G = max(0, (4 - abs(g(i) - 2) - abs(g(i) - 4))/2); 
    B = max(0, (3 - abs(g(i) - 1) - abs(g(i) - 2))/2); 
    cmapRainbow(i,:)=[R, G, B];
end

g=g';; %change g to column vector
f=f';; %change f to column vector


% % -----------------------------------------------
% % This program or any other program(s) supplied with it does not provide any warranty direct or implied.
% % This program is free to use/share for non-commercial purpose only. 
% % Kindly reference the work.
% % Author: Dr. Murtaza Khan
% % LinkedIn: http://www.linkedin.com/pub/dr-murtaza-khan/19/680/3b3
% % ResearchGate: http://www.researchgate.net/profile/Murtaza_Khan2/
% % -----------------------------------------------
