% t: current time, Y: y, p: params (type struct)
function dYdt = fn_immune_ode_system(t, Y, p, constMuE)
% Y order: [M0, M1, M2, K0, KA, E0, EF, EA]

    % Parameters --------------------------------------------
    [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
         KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma,rhoF, ...
            muE,omega,beta1,rho0] = default_parameters;
    
    % Input parameters
    omega = p.omega; % p("omega");
    beta1 = p.beta1; % p("beta1");
    rho0 = p.rho0; % p("rho0");

    % Calculate muE
    if (~constMuE)
        muE = (2*10^4).*(sin(pi*(t+12)/28)).^160;
    end

    % Calculate the system of DEs ----------------------------
    % Macrophage
    dM0dt = muM+etaM*(1-(Y(1)+Y(2)+Y(3))/MC)*(Y(1)+Y(2)+Y(3))-beta1*(Y(7)+Y(8))*Y(1)-thetaM*(Y(5)/(CKA+Y(5)))*Y(1)-beta2*Y(1)-deltaM*Y(1);
    dM1dt = beta1*(Y(7)+Y(8))*Y(1)+thetaM*(Y(5)/(CKA+Y(5)))*Y(1)+beta21*Y(3)-beta12*Y(2)-deltaM*Y(2);
    dM2dt = beta2*Y(1)+beta12*Y(2)-beta21*Y(3)-deltaM*Y(3);
    % Natural killer
    dK0dt = muK-thetaK*Y(4)*(Y(2)/(CM1+Y(2)))-deltaK*Y(4);
    dKAdt = thetaK*Y(4)*(Y(2)/(CM1+Y(2)))+etaK*(1-(Y(5)+Y(4))/KC)*Y(5)-sigma*(Y(7)+Y(8))*Y(5)-deltaK*Y(5);
    % Endometrial
    dE0dt = muE-deltaE0*(rho0*Y(6)+(1-rho0)*Y(6));
    dEFdt = deltaE0*rho0*Y(6)-omega*(gamma*Y(5)+(1-gamma)*Y(2))*Y(7)-rhoF*Y(7)-deltaE*Y(7);
    dEAdt = rhoF*Y(7)+etaE*(1-Y(8)/EC)*(Y(3)/(CM2+Y(3)))*Y(8)-omega*(gamma*Y(5)+(1-gamma)*Y(2))*Y(8)-deltaE*Y(8);

    % Concatenate into column array for return
    dYdt = [dM0dt; dM1dt; dM2dt; ...
            dK0dt; dKAdt; ...
            dE0dt; dEFdt; dEAdt];
end