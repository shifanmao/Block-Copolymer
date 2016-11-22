function f = fits_rm(x)
  global X1 Y1 X2 Y2 fitmin fitmax
  
  % define fitting parameters
  RM = x(1);
  SCALE = x(2);
  BG = x(3:end);

  % start parsing fits
  [~,NFIT] = size(X1);
  f = [];

  for IT = 1:NFIT
    x1 = X1(:,IT);y1 = Y1(:,IT);
    x2 = X2(:,IT);y2 = Y2(:,IT);

    % focus on fitting range
    xmin = fitmin; xmax = fitmax;
    minind = find(x1>=xmin);minind = minind(1);
    maxind = find(x1<=xmax);maxind = maxind(end);
    x1 = x1(minind:maxind);
    y1 = y1(minind:maxind);

    % target = log(y1);
    % Yfit = interp1(x2,log(y2),x1*RM);
    target = nearest_estimate(x1*RM,log(y1),x2,log(y2));
    func = SCALE*log(y2)+BG(IT);

    fdiff = target - func;
    fdiff(target==0) = 0;
    f = [f;fdiff];
  end
end

function y_est = nearest_estimate(x1,y1,x2,y2)
  y_est = zeros(length(x2),1);
  for ii = 1:length(y_est)
%    ind = find(abs(x1-x2(ii)) == min(abs(x1-x2(ii))));
%    y_est(ii) = y1(ind);

%    if ind==1 || ind==length(x1)
%      y_est(ii) = 0;
%    end
    y_est(ii) = interp1(x1,y1,x2(ii));
    if isnan(y_est(ii))
      y_est(ii) = 0;
    end
  end
end
