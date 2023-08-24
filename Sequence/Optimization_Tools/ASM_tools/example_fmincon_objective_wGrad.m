function [f,g] = example_fmincon_objective_wGrad(x)
% Calculate objective f
f = x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);

if nargout > 1 % gradient required
    grad_x1 = 2*x(1)*x(4) + x(2)*x(4) + x(3)*x(4);
    grad_x2 = x(1)*x(4);
    grad_x3 = x(1)*x(4) + 1;
    grad_x4 = x(1)^2 + x(1)*x(2) + x(1)*x(3);
    g       = [grad_x1; grad_x2; grad_x3; grad_x4];
end