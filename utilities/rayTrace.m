function inds = rayTrace(x1,x2)
% Trace a line through a 2d grid and find each pixel the line intersects
%
% http://playtechs.blogspot.com/2007/03/raytracing-on-grid.html

dx = abs(x1(1)-x2(1));
dy = abs(x1(2)-x2(2));
x = x1(1);
y = x1(2);
n = 1 + dx + dy;
if x1(1) > x2(1)
    x_inc = -1;
else
    x_inc = 1;
end

if x1(2) > x2(2)
    y_inc = -1;
else
    y_inc = 1;
end

error = dx - dy;
dx = 2*dx;
dy = 2*dy;

inds = zeros(n,2);
i = 1;
while i < n+1
    inds(i,:) = [x,y];
    
    if error > 0
        x = x + x_inc;
        error = error - dy;
    else
        y = y + y_inc;
        error = error + dx;
    end
    i = i + 1;
end

end