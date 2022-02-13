clear variables;
close all;

h = 0.01;
t = 0:h:2;
y0 = 1;

%calcul de y' = f(t, y)
%l'équa diff est (1+t*t)y' - 2ty = 3t-1
f = @(y,t) (2*t*y + 3*t -1)/(1+t*t);

y = zeros(4,length(t));
y(:,1) = y0;
beta = 0.5;

for k=1:length(t)-1
    %euler explicite
    y(1,k+1) = y(1,k) + h*f(y(1,k), t(k));
    
    %méthode d'euler implicite
    y(2,k+1) = ((1 + t(k+1)*t(k+1))*y(2,k) + h*(3*t(k+1) - 1))/(1 + t(k+1)*t(k+1) - 2*h*t(k+1));
    
    %runge-Kutta
    k1 = f(y(3,k),t(k));
    k2 = f(y(3,k) + (h*k1)/(2*beta)  ,t(k) + h/(2*beta));
    y(3,k+1) = y(3,k) + h*((1-beta)*k1 + beta*k2);
    
    %runge-kutta ordre 4
    k1 = f(y(4,k),t(k));
    k2 = f(y(4,k) + (h/2)*k1, t(k) + h/2);
    k3 = f(y(4,k) + (h/2)*k2, t(k) + h/2);
    k4 = f(y(4,k) + h*k3, t(k) + h);
    y(4,k+1) = y(4,k) + (h/6)*(k1 + 2*k2 +2*k3 + k4);
end

% affichage
figure(1);hold on;
plot(t,y(1,:),'c');
plot(t,y(2,:),'m');
plot(t,y(3,:),'r');
plot(t,y(4,:),'b');
lg=legend('Euler explicite','Euler implicite','RK2','RK4');