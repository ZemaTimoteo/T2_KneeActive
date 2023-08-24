clc, tic
% 1 - objective + gradients function
objective = @example_fmincon_objective_wGrad;

% 2 - initial guess
x0 = [1,5,5,1];

% 3 - variable bounds
lb = 1.0 * ones(4,1);
ub = 5.0 * ones(4,1);

% 4 - show initial objective
disp(['Initial Objective: ' num2str(objective(x0))])

% 5 - linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% 6 - nonlinear constraints
nonlincon = @example_fmincon_Matlab_aux;

% 7 - options
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);

% 8 - optimize with fmincon
x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon,options);

% 9 - show final objective
disp(['Final Objective: ' num2str(objective(x))])

% 10 - print solution
disp('Solution')
disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])
disp(['x4 = ' num2str(x(4))])
toc