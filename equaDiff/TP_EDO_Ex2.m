clear variables;
close all;

tmin=0;tmax=1;
f=@(t,y) t.^3*exp(-5*t) - (4*(t.^3) + 5)*y;
yExact=@(t) 0.25*(exp(t.^4) + 3).*exp(-t.*(t.^3+5));
    
% condition initiale
y0=1;

%% question 1
% 1. méthode d'Euler (h=0.1)
h=0.1;

[yEuler1,t1]=fct_Euler(y0,tmin,tmax,h,f);
eps1=abs(yEuler1-yExact(t1));   % erreur


% 2. méthode d'Euler (h=0.05)
h=0.05;
[yEuler2,t2]=fct_Euler(y0,tmin,tmax,h,f);
eps2=abs(yEuler2-yExact(t2));   % erreur

% 3. méthode RK2 (h=0.1 et beta=1)
h=0.1;beta=1;
[yRK,t3]=fct_RK2(y0,tmin,tmax,h,beta,f);
eps3=abs(yRK-yExact(t3));       % erreur

figure(1)
subplot(211)
plot(t1,yExact(t1),'r', t1, yEuler1, 'b', t2, yEuler2, 'g',t3, yRK)

subplot(212)
plot(t1,eps1, t2, eps2, t3,eps3 );


figure(2)

Thi = 0.02:0.001:0.1;

for k=1:length(Thi)
    %euler
    
    hi = Thi(k);
    [yEuler2,t2]=fct_Euler(y0,tmin,tmax,hi,f);
    eps2=abs(yEuler2-yExact(t2));   % erreur
    subplot(121)
    hold on;
    plot(hi,max(eps2),'*r')
    
    subplot(122)
    [yRK,t3]=fct_RK2(y0,tmin,tmax,hi,beta,f);
    eps3=abs(yRK-yExact(t3));       % erreur
    hold on;
    plot(hi,max(eps3),'*r')
    
end
