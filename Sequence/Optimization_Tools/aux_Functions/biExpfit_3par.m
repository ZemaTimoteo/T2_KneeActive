%==========================================================================
%   Bi exponential fitting
%   x: Time points.
%   y: array of signal intensities for each time point
%   y0: constant term to include noise
%   azadeh.sharafi@nyulabgone.org  last edited : 2016-09-21
%=========================================================================
function [shortT, longT,shortA,longA,fit_out,gof,resid]  = biExpfit_3par(x,y)

% 3 par fitting
g = fittype( 'a*exp(-x/b)+ (1-a)*exp(-x/c)','dependent','y', 'independent','x');
option = fitoptions(g);
option.Lower = [0 0 0]; %[0 0  0 0 ];
option.Upper = [1 40 200];%[];%
option.StartPoint = [0.5 5 50];
option.Algorithm = 'Trust-Region';%Levenberg-Marquardt';%
option.Robust = 'off';
option.Normalize = 'off';

%16.06.20
y = y./max(y(:));

[f,gof,o] = fit(x,y,g,option) ;
resid =o.residuals;
if f.b < f.c
    shortT = f.b;
    shortA = f.a;
    longA = 1-f.a;
    longT = f.c;
    
    %disp('Case f.b < f.c');
else
    shortT = f.c;
    shortA = 1-f.a;
    longA = f.a;
    longT = f.b;
    
     %disp('Case f.b > f.c');
end

% AF modified 30.03.2020
fit_out = f;
end