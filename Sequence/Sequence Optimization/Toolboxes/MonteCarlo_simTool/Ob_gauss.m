% Objective function to calculate difference F between Normal distribution,
% with mean and std given in vector x = [mu std], and a given signal y

function [F,jacF] = Ob_gauss(x,y,t)
% F = [x(1)*exp(-t/x(2))-y;0.01*x(1);0.01*x(2)];
% F = x(1)*exp(-t/x(2))-y;
F = 1/sqrt(2*pi)/x(2)*exp(-(t-x(1)).^2/(2*x(2)^2))-y;

if nargout > 1 % need Jacobian
    jacF = [exp(-t/x(2)) x(1)/x(2)^2.*t.*exp(-t/x(2));...
        0.01 0;0 0.01];
end