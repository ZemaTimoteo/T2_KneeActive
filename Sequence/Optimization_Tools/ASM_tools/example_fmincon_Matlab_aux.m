% create file nlcon.m for nonlinear constraints
function [c,ceq] = example_fmincon_Matlab_aux(x)
  c = 25.0 - x(1)*x(2)*x(3)*x(4);
  ceq = sum(x.^2) - 40;
