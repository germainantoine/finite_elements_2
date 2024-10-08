function val = g(x,y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
val=0;
if x==0 || x==2
    val = cos(2*pi*y);
end
if y==0 || y==2
    val = cos(pi*x);
end
end

