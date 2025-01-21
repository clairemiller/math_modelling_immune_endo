%% Function for determining the EA condition
function [lowEASystem, growth_rate, clear_rate] = ...
                                    fn_determine_condition_satisfied(T,parms)
        
    % Read in the default parameters
    [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
              KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, muE,omega,beta1,rho0] = default_parameters;
    
    % Overwrite the omega, beta1, and rho0 values
    omega = parms.omega;
    beta1 = parms.beta1;
    rho0 = parms.rho0;
    
    % Calculate the two sides of the equation
    M2_hat = T.M2./(CM2+T.M2);
    growth_rate = etaE*M2_hat;
    clear_rate = omega * (gamma*T.KA + (1-gamma)*T.M1) + deltaE;
    
    % Determine if the condition is true
    lowEASystem = (growth_rate < clear_rate);
end