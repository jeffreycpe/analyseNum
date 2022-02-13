clear variables;
close all;

h = 0.01;
t = 0:h:8;
y0 = 1;
y10 = 0;
y20 = 1;


e = @(t, y, u, v) u;
f = @(t, y, u, v) v;
g = @(t, y, u, v) t*exp(-t.^(2)) - 1/y ;

y = zeros(1, length(t));
y(1) = y0;
u = zeros(1, length(t));
u(1) = y10;
v = zeros(1, length(t));
v(1) = y20;
for k=1:length(t)-1
   y(k+1) = y(k) + h*e(t(k), y(k), u(k), v(k));
   u(k+1) = u(k) + h*f(t(k), y(k), u(k), v(k));
   v(k+1) = v(k) + h*g(t(k), y(k), u(k), v(k));
   
end 

figure(1)

plot(t,y)