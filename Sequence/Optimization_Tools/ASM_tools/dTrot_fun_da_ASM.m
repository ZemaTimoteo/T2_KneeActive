    function T = dTrot_fun_da_ASM(a,p)
        % derivative Rotation matrix w.r.t. alfa (tip angle amplitude)
        T = zeros([3 3]);
        T(1) = -sin(a)/2;%cos(a/2).^2;
        T(2) = 0.5*exp(-2*1i*p)*sin(a);%exp(-2*1i*p)*(sin(a/2)).^2;
        T(3) = -0.5*1i*exp(-1i*p)*cos(a);%-0.5*1i*exp(-1i*p)*sin(a);
        T(4) = conj(T(2));%conj(T(2));
        T(5) = T(1);%T(1);
        T(6) = 0.5*1i*exp(1i*p)*cos(a);%0.5*1i*exp(1i*p)*sin(a);
        T(7) = -1i*exp(1i*p)*cos(a);%-1i*exp(1i*p)*sin(a);
        T(8) = 1i*exp(-1i*p)*cos(a);%1i*exp(-1i*p)*sin(a);
        T(9) = -sin(a);%cos(a);
    end