function [] = circle(x,y,r)
    % Plot circle of center (x,y) and radius R
    %Input: x and y center
    %       r is radius
    %
    
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp, 'LineWidth', 2);

end