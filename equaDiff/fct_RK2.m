function [yRK,t]=fct_RK2(y0,tmin,tmax,h,beta,f)
    t=tmin:h:tmax;
    y = zeros(1,length(t));
    y(1) = y0;
    
    for k=1:length(t)-1
        k1 = f(t(k), y(k));
        k2 = f(t(k) + h/(2*beta), y(k) + (h*k1)/(2*beta));
        y(k+1) = y(k) + h*((1-beta)*k1 + beta*k2);
    end
    yRK = y;
end