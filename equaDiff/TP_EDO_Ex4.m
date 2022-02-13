clear variables;
close all;

tmin=0;
tmax=15;

% paramètres du modèle proies-prédateurs
alpha=1;  % taux de reproduction des proies
beta=0.5; % taux de mortalité des proies
gamma=2;  % taux de reproduction des prédateurs
delta=1;  % taux de mortalité des prédateurs

% conditions initiales
x0=2;     % proies
y0=0.5;   % predateurs

% second membre de l'équ. diff. (x',y')=(f(t,x,y),g(t,x,y))
f=@(t,x,y)(x.*(alpha-beta.*y));
g=@(t,x,y)(-y.*(gamma-delta.*x));


% Calcul des populations des proies et des prédateurs
h=0.01;   % pas temporel

% 1. méthode d'Euler
[xEuler,yEuler,t]=Euler_2D(x0,y0,tmin,tmax,h,f,g);



% 2. méthode RK2
beta=1;
[xRK2,yRK2,t]=RK2_2D(x0,y0,tmin,tmax,h,beta,f,g);


% Affichage des populations des proies et des prédateurs
% en fonction du temps
figure(1);hold on;

t = tmin:h:tmax;
subplot(121)
plot(t, xEuler, 'r', t, yEuler, 'b')
subplot(122)
plot(t, xRK2, 'r', t, yRK2, 'b')



% affichage de la trajectoire proies-prédateurs (tangente au champ de
% vecteurs défini par la fonction (x,y)->(f(t,x,y),g(t,x,y))
figure(2);hold on;

% champ de vecteurs (x,y)->(f(t,x,y),g(t,x,y))

N=40;
ux=linspace(0,7,N);
uy=linspace(0,7,N);
[x,y]=meshgrid(ux,uy);       % grille de coordonnées (ux,uy)
fxy=f(t,x,y);gxy=g(t,x,y);   % calcul du champ de vecteurs
norme=(fxy.^2+gxy.^2).^0.5;  % normalisation des vecteurs
fxy=fxy./norme;gxy=gxy./norme;
quiver(x,y,fxy,gxy);         % affichage de fxy et gxy sous forme
%                              % de champ de vecteurs


% trajectoires proies-prédateurs
% 1. méthode de d'Euler
plot(xEuler, yEuler, 'r')


% 2a. méthode RK2 
figure(3);hold on;
quiver(x,y,fxy,gxy);
plot(xRK2, yRK2, 'r')

% 2b. méthode RK2 et nouvelles conditions initiales
x0 = 3; 
y0 = 2;
[xRK2a,yRK2a,t] = RK2_2D(x0,y0,tmin,tmax,h,beta,f,g);
hold on;
plot(xRK2a,yRK2a, 'g')


