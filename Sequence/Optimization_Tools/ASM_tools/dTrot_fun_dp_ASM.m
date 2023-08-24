    function T = dTrot_fun_dp_ASM(a,p)
        % derivative Rotation matrix w.r.t. p (tip angle phase)
        T = zeros([3 3]);
        T(1) = 0;%cos(a/2).^2;
        T(2) = -2*1i*exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));
        T(5) = T(1);%T(1);
        T(6) = -0.5*exp(1i*p)*sin(a);
        T(7) = exp(1i*p)*sin(a);
        T(8) = exp(-1i*p)*sin(a);
        T(9) = 0;%cos(a);
    end