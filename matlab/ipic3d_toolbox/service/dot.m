function dot = dot(ax, ay, az, bx, by, bz, n)
dot = ax.*bx + ay.*by + az.*bz;
if (nargin>6)
dot = dot ./n;
end 
