clc, tic
% 1 - objective function
objective = @(x) x(1)*x(4)*(x(1)+x(2)+x(3))+x(3);

% 2 - initial guess
x0 = [1,5,5,1]; 
b = [];
Aeq = [];
beq = [];

% 6 - nonlinear constraints
nonlincon = @example_fmincon_Matlab_aux;

% 7 - optimize with fmincon
%[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]
% = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon); % x - VFA's

% 8 - show final objective
disp(['Final Objective: ' num2str(objective(x))])

% 9 - print solution
disp('Solution')
disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])
disp(['x4 = ' num2str(x(4))])
toc