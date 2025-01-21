function [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
              KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, muE,omega,beta1,rho0] = default_parameters
    muM = 1e4;
    muK = 3.5e3; % old: 1e4, new: 3.5e3
    etaE = 1;
    etaK = 0.02;
    etaM = 0.2;
    deltaM = 0.02;
    deltaK = 0.02;
    deltaE0 = 1;
    deltaE = 0.14;
    sigma = 1e-5;
    MC = 1e6;
    KC = 1e6;
    EC = 1e9;
    CM1 = 1e5;
    CM2 = 100;
    CKA = 1e5; % old: 1e6, new: 1e5
    beta2 = 1e-3;
    beta12 = 5e-5;
    beta21 = 5e-5;
    thetaM = 5e-8;
    thetaK = 0.7; % old: 0.2, new: 0.7;
    gamma = 0.8;
    rhoF = 0.1;
    muE = 1260;
    omega = 1e-5;
    beta1 = 1e-6;
    rho0 = 0.1;
end
