function [yEuler,t]=fct_Euler(y0,tmin,tmax,h,f)
    
    t=tmin:h:tmax;
    y = zeros(1,length(t));
    y(1) = y0;
    
    for k=1:length(t)-1
        y(k+1) = y(k) + h*f(t(k), y(k));
    end
    yEuler = y;
end