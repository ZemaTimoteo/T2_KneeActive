function [VarT2, epg] = Constraint_test_CRLB(x,ETL,params)
% Cost Function and Gradient obtained for the variance expected for a specific model/set of parameters
%
% Functions used:
%   Constraint_PS_CRLB: to generate the constraint of CRLB for
%   Patternsearch
%
% Inputs:
%   x = [TE (ms), FA(rad)]
%   ETL
%   params - parameters
% Ouputs:
%   Constraint of CRLB: Variance with respect to T2 estimation - Cannot be
%       smaller than
% 
%  @TiagoTFernandes, IST, Oct23

%% ... 0 - Parameters ...

% 0.1 - Inputs
TE    = x(1);              % beta = TE - Echo Time (ms)
FA    = x(2:end);          % Flip Angle constant (rad) 

% 0.2 - Other parameters
res    = params.res * params.accFactor;     % Resolution by Accelerator Factor (only acquired the ones needed - GRAPPA)
nsli   = params.nsli;                       % Number of slices
B1     = params.B1;                         % B1 value - TODO change
T1     = params.T1;                         % T1 (ms)
T2     = params.T2;                         % T2 (ms)
SNR    = params.SNR;     
sigma1 = params.sigma1;                     % RFexcDur in (ms) - taken from pulseq implementation of a 90º RF exc pulse
sigma2 = params.sigma2;                     % t_gs4 + t_gs5 + t_spoiler in (ms) - time of spoiler gradients
sigma3 = params.sigma3;                     % sigma of SNR

% 0.3 - Parameters for derivatives - Initialize
beta     = TE;              % interecho spacing in (ms) - Full of echospacing (ESP) in (ms) = TE
alpha_RF = params.alpha_RF; % Flip Angle in (rad)


%% 1 - Get gamma and d_gamma
% 1.1 - gamma thetformulation
numer_gamma_beta = (  1 - exp( - ( ( sigma1/2 + (1/2+ETL)*beta )*(nsli-1)+sigma2*nsli) / T1 ) )  ^ 2;
denom_gamma_beta = (sigma3^2/res) * sqrt( (sigma1/2 + ( 1/2 + ETL )*beta + sigma2)*nsli*res);
gamma.gamma_beta = numer_gamma_beta / denom_gamma_beta;

% 1.2 - Derivative of gamma in order to beta (ms)
d_numer_gamma_beta_dbeta = (  2*(1/2 + ETL)    *    (-1 + nsli)  *  exp(  (  -nsli*sigma2 - (-1 + nsli)*( sigma1/2 + (1/2 + ETL)*beta ) ) / T1  )  * ...
                              ( 1 - exp(  (-nsli*sigma2 - (-1 + nsli)*(sigma1/2 + (1/2 + ETL)*beta))/T1  ) ) ...
                             ) /  T1;
d_denom_gamma_beta_dbeta = ((1/2 + ETL) * nsli * sigma3^2) / (  sqrt(2)*sqrt(nsli*res * (sigma1 + 2*sigma2 + beta + 2*ETL*beta))  );


gamma.dgamma_beta = ( denom_gamma_beta * d_numer_gamma_beta_dbeta  - numer_gamma_beta * d_denom_gamma_beta_dbeta) / ...
                         (denom_gamma_beta^2)  ;

                     
%% ... 2 - Calculate Gradients (derivatives) of EPG
% 2.1 - Dictionary with SLR Profile
if params.methodDic == 'SLR_Prof'    
    clear grad_tool FF_tool
    % 2.1.1 - RF pulses
    refoc_puls = [];
    for jj=1:size(FA,2)
        exc_puls                   = params.alpha_exc; % Flip Angle in (rad)       
        degree_FA                  = FA(jj)*180/pi;
        [exc_puls, aux_refoc_puls] = slr_profile_wRFoptimz_design(exc_puls,B1,degree_FA,params);
        refoc_puls                 = [refoc_puls aux_refoc_puls(:)];  
    end

    % 2.1.2 - Get Derivatives
    s_epg = [];
    s_grad_ds_dT2 = [];
    for z=1:size(refoc_puls,1)
        [aux_grad, aux_x ] = CF_Mycode_epg_derivatives_NLO_AUC(...
                                ETL, beta, ...
                                T1, T2, exc_puls(z), ...
                                refoc_puls(z,:), params, gamma); % [dF_dalphaL,dF_dbeta]
        s_epg(:,z)         = aux_x;
        ds_dT2(:,z)         = aux_grad.ds_dT2;
    end
    
    grad.ds_dT2           = sum(ds_dT2,2);
    
% 2.2 - Just Dictionary                        
elseif params.methodDic == 'JUSTdict' 
    exc_puls                   = params.alpha_exc; % Flip Angle in (rad)
    [grad, x ]= CF_Mycode_epg_derivatives_NLO_AUC(...
                            ETL, beta, ...
                            T1, T2, exc_puls, ...
                            FA, params, gamma); % [dF_dalphaL,dF_dbeta]  
end


ds_dT2 = grad.ds_dT2; % derivative over T2


%% ... 3 - CRLB & Signal ...
    % --- 3.1 Calculate CF value - Uncertainty of CRLB ---
data.vardT2 = gamma.gamma_beta * (ds_dT2(:)'*ds_dT2(:));

%% ... 4 - CRLB Constraint
% Cannot be samller then (so it has to be negative cause result as to be bigger then standart)
%                | maxVart2 with respect to specific T2 in constant FlipAngle approach
VarT2 = data.vardT2;

end

