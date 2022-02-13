function [x,y,t]=RK4_2D(x0,y0,tmin,tmax,pas,beta,F,G)
    h=pas;
    t=tmin:h:tmax;
    y = zeros(1,length(t));
    x = zeros(1,length(t));
    x(1) = x0;
    y(1) = y0;
    beta = 1;
    for k=1:length(t)-1
        k1x = F(t(k), x(k), y(k));
        k1y = G(t(k), x(k), y(k));
        
        k2x = F(t(k) + h/(2*beta), x(k) + (h*k1x)/(2*beta), y(k) + (h*k1y)/(2*beta));
        k2y = G(t(k) + h/(2*beta),x(k) + (h*k1x)/(2*beta), y(k) + (h*k1y)/(2*beta));
        
        k3x = F(t(k) + h/(2*beta), x(k) + (h*k2x)/(2*beta), y(k) + (h*k2y)/(2*beta));
        k3y = G(t(k) + h/(2*beta),x(k) + (h*k2x)/(2*beta), y(k) + (h*k2y)/(2*beta));
        
        k4x = F(t(k) + h, x(k) + (h*k3x), y(k) + (h*k3y));
        k4y = G(t(k) + h, x(k) + (h*k3x), y(k) + (h*k3y));
        
        x(k+1) = x(k) + (h/6)*(k1x + 2*k2x + 2*k3x + k4x);
        y(k+1) = y(k) + (h/6)*(k1y + 2*k2y + 2*k3y + k4y);
    end

  end