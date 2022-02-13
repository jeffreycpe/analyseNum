clear variables;
close all;

% paramètres physiques
m=0.7; % masse de la bille (kg)
r=0.035; % rayon de la bille (m)
eta=0.000018; % coeff. de viscosité de l'air à 20°C (kg.m^-1.s^-1)
gamma=6*pi*r*eta/m; % frottements (s^-1)
gr=9.8; % accéleration de la pesanteur (m.s^-2)
l=2; % longueur du fil (m)
omega=sqrt(gr/l); % fréquence propre (rad.s^-1)
T0=2*pi/omega; % (pseudo-)période du pendule (s)

% autres paramètres
tmin=0; % instant initial
tmax=4*T0; % instant final
pas=0.001; % pas de calcul
%fprintf('Durée de l''expérience physique : %1.2f\n',tmax-tmin);


% fonctions Y'=F(Y) avec ici Y=(theta,z) et F(Y)=(f,g)
f=@(t,teta,z)(z);
g=@(t,teta,z)(-gamma*1000*z - (omega^2)*sin(teta));


% conditions initiales
theta0=2*pi/3; % angle initial (rad)
thetap0=0; % vitesse angulaire initiale (rad/s)

%euler
[theta,thetap,t]=Euler_2D(theta0,thetap0,tmin,tmax,pas,f,g);

% 2. méthode RK2
beta=1;
[RKtheta,RKthetap,RKt]=RK2_2D(theta0,thetap0,tmin,tmax,pas,beta,f,g);

%méthode RK24

beta=1;
[RK4theta,RK4thetap,RK4t]=RK4_2D(theta0,thetap0,tmin,tmax,pas,beta,f,g);

figure(1)
subplot(221)
plot(theta, thetap);
subplot(222)
plot(RKtheta, RKthetap);
subplot(223)
plot(RK4theta, RK4thetap);

figure(2)
Ec = 0.5*m*l*l*(thetap.^2);
Ep = m*gr*l*(1 - cos(theta));
Em = Ec + Ep;

subplot(221)
plot(t, Ec, t, Ep, t, Em,'r')

subplot(222)
RKEc = 0.5*m*l*l*(RKthetap.^2);
RKEp = m*gr*l*(1 - cos(RKtheta));
RKEm = RKEc + RKEp;

plot(t, RKEc, t, RKEp, t, RKEm,'r')

subplot(223)
RK4Ec = 0.5*m*l*l*(RK4thetap.^2);
RK4Ep = m*gr*l*(1 - cos(RK4theta));
RK4Em = RK4Ec + RK4Ep;

plot(t, RK4Ec, t, RK4Ep, t, RK4Em,'r')