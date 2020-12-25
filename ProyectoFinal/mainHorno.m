clc; clear; close all;

%% Tiempo de simulacion
time = [0:0.001:10];

%% parametros planta horno
params.V = 10;
params.Tin = 298;
params.Top = 298;
params.rho_a = 1;
params.UA = 1;
params.rho_a = 1;
params.Cp = 1;

%% Condiciones iniciales
% X0 = [CA;CB;CC;CD;CE;CF;T;Fin];
X0 = [1;2;3;4;5;6;7;8];

%% Simulacion
% horno(X, params)
[t, result] = ode45(@horno, time, X0, [], params);

figure();
for i = 1:8
    subplot(2,4,i)
    plot(t, result(:,i),'LineWidth',2);
    xlabel('[s]')
    grid on
end

function Y = horno(time, X, params)
    var = num2cell(X);
    [Ca, Cb, Cc, Cd, Ce, Cf, T, Fin] = var{:};
    
    V = params.V;
    Tin = params.Tin;
    Top = params.Top;
    UA = params.UA;
    rho_a = params.rho_a;
    Cp = params.Cp;
    
    dT_dt = (Fin*(Tin - T)/V) - (UA*(Top - T)/(rho_a*V*Cp));
    Y = [0; 0; 0; 0; 0; 0; dT_dt; 0];
end