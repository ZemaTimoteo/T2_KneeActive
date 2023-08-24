%==========================================================================
%   Mono exponential fitting
%   x: Time points.
%   y: array of signal intensities for each time point
%   y0: constant term to include noise
%   azadeh.sharafi@nyulabgone.org  last edited : 2016-09-21
%=========================================================================
function [singleT,fit_out,A,gof,resid]  = monoExpfit(x,y)

g = fittype( 'a*exp(-x/b)', 'independent','x');
option = fitoptions(g);
option.Lower = [0 0]; %[0 0  0 0 ];
option.Upper = [Inf 500];%[];%
option.StartPoint = [1 5];
option.Algorithm = 'Trust-Region';%Levenberg-Marquardt';%
option.Robust = 'off';
option.Normalize = 'off';

%16.06.20
y = y./max(y(:));

[f,gof,o] = fit(x,y,g,option) ;
resid =o.residuals;
singleT = f.b;
A = f.a;


% AF modified 23.03.2020
fit_out = f;

end
        
        