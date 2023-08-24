%% Calculates CRLB value for EPG of FSE

function f = CRLB_EPG(x,ETL,T2)

    t = x:x:ETL*x;
    M0 = 1;
    sigma= 1e-2;

    % S = M0*exp(-t./T2);

    dS_dM0 = exp(-t./T2);
    dS_dT2 = M0.*(t./T2^2).*exp(-t./T2);

    dM      = [dS_dM0(:) dS_dT2(:)];

    FIM    = dM_num.'   *  (   (1/sigma^2) * eye( size(dM_num,1) )  )   *  dM_num;
    CRLB   = diag(  sqrt( inv(FIM) )  );


    f = CRLB(2)./T2;  % normalise the CRLB of each parameter by its respective value.


end