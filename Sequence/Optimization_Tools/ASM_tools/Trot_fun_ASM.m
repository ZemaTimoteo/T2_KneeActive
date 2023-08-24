% Rotation matrix direct definition
function T = Trot_fun_ASM(a,p)

T = zeros([3 3]);
T(1) = cos(a/2).^2;
T(2) = exp(-2*1i*p)*(sin(a/2)).^2;
T(3) = -0.5*1i*exp(-1i*p)*sin(a);
T(4) = conj(T(2));
T(5) = T(1);
T(6) = 0.5*1i*exp(1i*p)*sin(a);
T(7) = -1i*exp(1i*p)*sin(a);
T(8) = 1i*exp(-1i*p)*sin(a);
T(9) = cos(a);
end