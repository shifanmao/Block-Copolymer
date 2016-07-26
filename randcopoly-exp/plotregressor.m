function [xp,yp] = plotregressor(x,y,fitrange,plotrange)
  xfit = fitrange;  % linear fit regime
  xf = x(xfit)';
  Xp = [ones(length(xf),1) xf];
  yp = y(xfit)';
  b = Xp\yp;

  xp = linspace(x(plotrange(1)),x(plotrange(end)),10);
  yp = b(2)*xp+b(1);
end
