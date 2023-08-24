function E = dErelax_fun_dT2_ASM(T2,TE,S0,N)
    % derivative Relaxation matrix w.r.t. T2
    E  = eye(N);
    dE1 = 0;                              % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF
    dE2 = S0*((TE)/(T2.^2)).*exp(-TE/T2); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF
% %     dE2 = S0 * ( (-TE) ) * exp(-TE/T2); % dividing ESP by to i.o.t. be similar of epg_cpmg_ASM - to split before and after RF

    ii                       = 1:(N/3);
    E(3*N*(ii-1)+3*(ii-1)+1) = dE2;
    E(3*N*ii-2*N+3*(ii-1)+2) = dE2;
    E(3*N*ii-N+3*(ii-1)+3)   = dE1;
end