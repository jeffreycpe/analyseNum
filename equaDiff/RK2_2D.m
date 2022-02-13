function [x,y,t]=RK2_2D(x0,y0,tmin,tmax,pas,beta,F,G)
    h=pas;
    t=tmin:h:tmax;
    y = zeros(1,length(t));
    x = zeros(1,length(t));
    x(1) = x0;
    y(1) = y0;
    
    for k=1:length(t)-1
        k1x = F(t(k), x(k), y(k));
        k1y = G(t(k), x(k), y(k));
        
        k2x = F(t(k) + h/(2*beta), x(k) + (h*k1x)/(2*beta), y(k) + (h*k1y)/(2*beta));
        k2y = G(t(k) + h/(2*beta),x(k) + (h*k1x)/(2*beta), y(k) + (h*k1y)/(2*beta));
        
        x(k+1) = x(k) + h*((1-beta)*k1x + beta*k2x);
        y(k+1) = y(k) + h*((1-beta)*k1y + beta*k2y);
    end

  end