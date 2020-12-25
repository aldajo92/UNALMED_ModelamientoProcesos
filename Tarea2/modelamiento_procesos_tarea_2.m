clc; clear; close all

% Condiciones iniciales
X0 = [1.305;0;0;0;0;350;420];

plant(0, X0)

time = [0 100];
[t, result] = ode45(@plant, time, X0, []);

figure()

n = 7;
m = 1;

subplot(n,m,1)
plot(t, result(:,1),'LineWidth',2);
title('CA')
xlabel('time (s)')
ylabel('[Kmol/m^3]')

subplot(n,m,2)
plot(t, result(:,2),'LineWidth',2);
title('CT')
xlabel('time (s)')
ylabel('[Kmol/m^3]')

subplot(n,m,3)
plot(t, result(:,3),'LineWidth',2);
title('CR')
xlabel('time (s)')
ylabel('[Kmol/m^3]')

subplot(n,m,4)
plot(t, result(:,4),'LineWidth',2);
title('CS')
xlabel('time (s)')
ylabel('[Kmol/m^3]')

subplot(n,m,5)
plot(t, result(:,5),'LineWidth',2);
title('CU')
xlabel('time (s)')
ylabel('[Kmol/m^3]')

subplot(n,m,6)
plot(t, result(:,6),'LineWidth',2);
title('T')
xlabel('time (s)')
ylabel('[K]')

subplot(n,m,7)
plot(t, result(:,7),'LineWidth',2);
title('Tj')
xlabel('time (s)')
ylabel('[K]')

function Y = plant(t, X)
% X : [CA,CT,CR,CS,CU,T,Tj]

    CA = X(1);
    CT = X(2);
    CR = X(3);
    CS = X(4);
    CU = X(5);
    T = X(6);
    Tj = X(7);
    
    if (CA < 0)
        CA = 0;
    end
    
    if (CT < 0)
        CT = 0;
    end
    
    if (CR < 0)
        CR = 0;
    end
    
    if (CS < 0)
        CS = 0;
    end
    
    if (CU < 0)
        CU = 0;
    end
    
    Runiv = 8.314472;   % [J/(K mol)] : Constante universal de los gases
    
    A0_vec = [
        5.8e11;         % [litro/(Kmol·s)]
        1.6e10;         % [s^(-1)]
        3.6e3;          % [s^(-1)]
        4.8e8           % [litro/(Kmol·s)]
    ];
    Ea_vec = [
        46275;          % [J/mol]
        59412;          % [J/mol]
        23137;          % [J/mol]
        30480           % [J/mol]
    ];
    delta_H = [
        -31.79;         % [KJ/kmol]
        2.45;           % [KJ/kmol]
        -25.56;         % [KJ/kmol]
        1.38            % [KJ/kmol]
    ];
    
    k = A0_vec .* exp(-Ea_vec/(Runiv*T));
    
    g = 9.8;            % [m/s^2]
    
    wA = 1.305*g;       % [kg/(m^2 s^2)]
    wT = 1.5*g;         % [kg/(m^2 s^2)]
    wR = 1.458*g;       % [kg/(m^2 s^2)]
    wS = 1.45*g;        % [kg/(m^2 s^2)]
    wU = 1.83*g;        % [kg/(m^2 s^2)]
    
    CA0 = 1.305;        % [Kmol/m^3]
    F = 140;            % [m^3/s]
    V = 1570;           % [m^3]
    Vchaqueta = 150;    % [m^3]
    Cp = 11.6162;       % [Kj/kg ºK]
    Cpj = 4.182;        % [Kj/kg ºK]
    ro = 876;           % [kg/m^3] : densidad solvente (benceno)
    roH2O = 1000;       % [kg/m^3] : densidad refrigerante (agua)
    U = 23;             % [Kj/m^2 ºK] coeficiente global de transferencia
    A = 236;            % [m^2]
    
    TA = 350;           % [K]
    Tj0 = 420;          % [K]
    
    dCA_dt = ((CA0*F)/V) - ((CA*F)/V) - (k(1)*(CA^2)/wA) - (k(2)*CA);
    dCT_dt = (wT*k(2)*(CA)/wA) - (CT*F/V);
    dCR_dt = (wR*(k(1)*(CA/wA)^2)/2) - (k(3)*(CR)) - (k(4)*(CR^2)/wR) - (CR*F/V);
    dCS_dt = (wS*k(3)*CR/wR) - (CS*F/V);
    dCU_dt = (wU*k(4)*((CR/wR)^2)/2) - (CU*F/V);
    
    r = [k(1)*(CA/wA)^2; k(2)*CA/wA; k(3)*CR/wR; k(4)*(CR/wR)^2];
    
    Q = U*A*(T-Tj);
    
    dT_dt = F*((CA0*TA) - CA*T)/(ro*V) + ((r'*delta_H)/(ro*Cp)) - (Q/(V*ro*Cp));
    dTj_dt = (F/V)*(Tj0 - Tj) + (Q/(Vchaqueta*roH2O*Cpj));
    
    Y = [dCA_dt; dCT_dt; dCR_dt; dCS_dt; dCU_dt; dT_dt; dTj_dt];
end
