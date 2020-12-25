clc; clear; close all
format SHORTG
%% Tiempo de simulacion
time = [0 50];

syms Ca Cb Cc Cd Ce T1 Fin1 Fin2
% Fin1 = 1176;
% Fin2 = 185.5;
%% Parametros de la planta
V = 2000;   % [] Volumen del reactor
Vc = 500;  % [] Volumen del condensador
R = 0.082;  % [atm lt/mol K] Constante de los gases ideales
UA = 34;    % [KJ/Ks]
UAcond = 20;    % [KJ/Ks]
PH_2 = 48;  % [atm] Presión del hidrogeno en la entrada
TH_2 = 298; % [K] Temperatura del hidrógeno en la entrada
wH_2 = 2;   % [gr/mol] peso molecular del hidrogeno
ro_mezcla = 2.35; % [] Densidad gases entrada
ro_h2 = PH_2*wH_2/(R*TH_2);  % [gr / lt] Densidad del hidrogeno
ro_B = PH_2*wH_2/(R*T1);  % [gr / lt] Densidad del hidrogeno
Ta = 295;       % [K] Temperatura ambiente
Tgases0 = 400;  % [K] Temperatura entrada Gases (mezcla)
Tcond = 337.8;  % [K] Temperatura condensador
Tref = Ta;

Cpgases = 0.0035;   % [KJ/mol K] Calor especifico gases
CpA = 0.03712;      % [KJ/mol K] Calor especifico monoxido de carbono
CpB = 0.0288;       % [KJ/mol K] Calor especifico hidrogeno
CpC = 0.081;        % [KJ/mol K] Calor especifico metanol
CpD = 0.01922;      % [KJ/mol K] Calor especifico dioxido de carbono
CpE = 0.0745;       % [KJ/mol K] Calor especifico agua

deltaHr1 = 128.2;   % [KJ/mol]
deltaHr2 = -41.2;   % [KJ/mol]
H1 = CpB * (Ta - Tref);
H2 = Cpgases * (Tgases0-Tref);

Ha = CpA * (T1 - Ta);
Hb = CpB * (T1 - Ta);
Hc = CpC * (T1 - Ta);
Hd = CpD * (T1 - Ta);
He = CpE * (T1 - Ta);

%Fracciones masicas
ya = 0.6; % Fraccion masica CO
yd = 0.3; % Fraccion masica CO2

wa = 28.01;
wb = 2;
wd = 44.01;

%Condiciones de entrada
Ca0 = ro_mezcla * ya/wa;
Cb0 = ro_h2/wb;
Cd0 = ro_mezcla * yd/wd;

Fout = Fin1 + Fin2;

k1 = 0.933*exp(2.5*(31400/1.987)*((1/330) - (1/T1)));
k2 = 0.933*exp((18000/1.987)*((1/300) - (1/T1)));

r_1 = k1*Ca*(Cb^2);
r_2 = k2*Cd*Cb;

r_c = r_1;
r_a = r_2;
r_e = r_2;

r_atoc = r_c;
r_btoc = 2*r_c;
r_btoa = r_a;
r_dtoa = r_a;

%% Ecn dinamicas del reactor
sumDeltaHrs = V*(r_1*deltaHr1 + r_2*deltaHr2);
Q = 34*(300-Tref);
sumCiHi = Ca*Ha + Cb*Hb + Cc*Hc + Cd*Hd + Ce*He;

dCa_dt = (Fin1*Ca0/V) - (Fout*Ca/V) - r_atoc + r_a;
dCb_dt = (Fin2*Cb0/V) - (Fout*Cb/V) - r_btoc + r_btoa;
dCc_dt = -(Fout*Cc/V) + r_c;
dCd_dt = (Fin1*Cd0/V) - (Fout*Cd/V) - r_dtoa;
dCe_dt = -(Fout*Ce/V) + r_e;
dT1_dt_num = (Fin1*ro_B*H1) + (Fin2*ro_mezcla*H2) - (Fout*sumCiHi) + Q - sumDeltaHrs;
dT1_dt_den = (ro_mezcla*Cpgases + ro_B*CpB)*V;
dT1_dt = dT1_dt_num/dT1_dt_den;

F = [dCa_dt; dCb_dt; dCc_dt; dCd_dt; dCe_dt; dT1_dt];
%% Definicion de variables de estado, entradas y salidas
X = [Ca; Cb; Cc; Cd; Ce; T1];
U = [Fin1; Fin2];
Y = [Fout*Cc];
%% Se emplea el solver para obtener el punto de equilibrio
res = solve(F, X);
op = [res.Ca res.Cb res.Cc res.Cd res.Ce res.T1];
% X0_1 = double(op');
%% punto de equilibrio
X0_1 = [0.05225 0.27364 0.00277 0.00230 0.01153 328.18469];
plant(0, X0_1)
%% Evaluacion punto de equilibrio
plant(0, X0_1)
%% Jacobianos para la linealización
Al = jacobian(F,X);
Bl = jacobian(F,U);
Cl = jacobian(Y,X);
Dl = jacobian(Y,U);

%% Evaluamos en el punto de equilibrio los jacobianos
Ca = X0_1(1);
Cb = X0_1(2);
Cc = X0_1(3);
Cd = X0_1(4);
Ce = X0_1(5);
T1 = X0_1(6);

Fin10 = 1176;
Fin20 = 185.5;
Fout0 = Fin10 + Fin20;

Fin1 = Fin10;
Fin2 = Fin20;

A = eval(Al)
B = eval(Bl)
C = eval(Cl)
D = eval(Dl)

%% Analisis de estabilidad
eig_valores_A = eig(A)

% X0_1 = [Cs0;Cxeq1]; % result = [64; 0]

% Fin1 = 1176;
% Fin2 = 185.5;
% 
% X0_1 = [0 0 0 0 0 0];
% [t, result] = ode45(@plant, time, X0_1, []);
% 
% figure()
% for i = 1:6
%     subplot(2,3,i)
%     plot(t, result(:,i),'LineWidth',2);
%     xlabel('[h]')
%     ylabel('[gr/lit]')
%     grid on
% end
% 
% plant(0, [0.050339 1 0 0.016019 0 298 0 0]);



X0_1 = [0.050339 3.9286 0 0.016019 0 337.8];
[t, result] = ode45(@plant, time, X0_1, []);

figure()
variable = arrayfun(@char, X, 'uniform', 0);
for i = 1:6
    subplot(2,3,i)
    plot(t, result(:,i),'LineWidth',2);
    xlabel('[s]')
    title(string(variable(i)))
    grid on
end

figure()
out = Fout0*result(:,3);
plot(t, out,'LineWidth',2)
xlabel('[s]')
title("Productividad C_{CH_3OH}")
grid on


%% planta
function Y = plant(t, X)
% X : [Ca;Cb;Cc;Cd;T1;T2]

    Ca = X(1);
    Cb = X(2);
    Cc = X(3);
    Cd = X(4);
    Ce = X(5);
    T1 = X(6);
%     yc = X(7);
%     ye = X(8);
    Fin1 = 1176;
    Fin2 = 185.5;
    
    %Parametros de la planta
    V = 2000;   % [] Volumen del reactor
    Vc = 500;  % [] Volumen del condensador
    R = 0.082;  % [atm lt/mol K] Constante de los gases ideales
    UA = 34;    % [KJ/Ks]
    UAcond = 20;    % [KJ/Ks]
    PH_2 = 48;  % [atm] Presión del hidrogeno en la entrada
    TH_2 = 298; % [K] Temperatura del hidrógeno en la entrada
    wH_2 = 2;   % [gr/mol] peso molecular del hidrogeno
    ro_mezcla = 2.35; % [] Densidad gases entrada
    ro_h2 = PH_2*wH_2/(R*TH_2);  % [gr / lt] Densidad del hidrogeno
    ro_B = PH_2*wH_2/(R*T1);  % [gr / lt] Densidad del hidrogeno
    Ta = 295;       % [K] Temperatura ambiente
    Tgases0 = 400;  % [K] Temperatura entrada Gases (mezcla)
    Tcond = 337.8;  % [K] Temperatura condensador
    Tref = Ta;
    
    Cpgases = 0.0035;   % [KJ/mol K] Calor especifico gases
    CpA = 0.03712;      % [KJ/mol K] Calor especifico monoxido de carbono
    CpB = 0.0288;       % [KJ/mol K] Calor especifico hidrogeno
    CpC = 0.081;        % [KJ/mol K] Calor especifico metanol
    CpD = 0.01922;      % [KJ/mol K] Calor especifico dioxido de carbono
    CpE = 0.0745;       % [KJ/mol K] Calor especifico agua
    
    deltaHr1 = 128.2;   % [KJ/mol]
    deltaHr2 = -41.2;   % [KJ/mol]
    H1 = CpB * (Ta - Tref);
    H2 = Cpgases * (Tgases0-Tref);
    
    Ha = CpA * (T1 - Ta);
    Hb = CpB * (T1 - Ta);
    Hc = CpC * (T1 - Ta);
    Hd = CpD * (T1 - Ta);
    He = CpE * (T1 - Ta);
    
    %Fracciones masicas
    ya = 0.6; % Fraccion masica CO
    yd = 0.3; % Fraccion masica CO2
    
    wa = 28.01;
    wb = 2;
    wd = 44.01;
    
    %Condiciones de entrada
    Ca0 = ro_mezcla * ya/wa;
    Cb0 = ro_h2/wb;
    Cd0 = ro_mezcla * yd/wd;
    
    Fout = Fin1 + Fin2;
    
    k1 = 0.933*exp(2.5*(31400/1.987)*((1/330) - (1/T1)));
    k2 = 0.933*exp((18000/1.987)*((1/300) - (1/T1)));
    
    r_1 = k1*Ca*(Cb^2);
    r_2 = k2*Cd*Cb;
    
    r_c = r_1;
    r_a = r_2;
    r_e = r_2;
    
    r_atoc = r_c;
    r_btoc = 2*r_c;
    r_btoa = r_a;
    r_dtoa = r_a;
    
    % Para el reactor
    sumDeltaHrs = V*(r_1*deltaHr1 + r_2*deltaHr2);
    Q = 34*(300-Tref);
    sumCiHi = Ca*Ha + Cb*Hb + Cc*Hc + Cd*Hd + Ce*He;
    
    dCa_dt = (Fin1*Ca0/V) - (Fout*Ca/V) - r_atoc + r_a;
    dCb_dt = (Fin2*Cb0/V) - (Fout*Cb/V) - r_btoc + r_btoa;
    dCc_dt = -(Fout*Cc/V) + r_c;
    dCd_dt = (Fin1*Cd0/V) - (Fout*Cd/V) - r_dtoa;
    dCe_dt = -(Fout*Ce/V) + r_e;
    dT1_dt_num = (Fin1*ro_B*H1) + (Fin2*ro_mezcla*H2) - (Fout*sumCiHi) + Q - sumDeltaHrs;
    dT1_dt_den = (ro_mezcla*Cpgases + ro_B*CpB)*V;
    dT1_dt = dT1_dt_num/dT1_dt_den;
    
    % Para el condensador
%     NC = Vc*Cc
%     NE = Vc*Ce
%     NCE = NC+NE
%     n_C_in = Fout * Cc;
%     n_E_in = Fout * Ce;
%     n_CE_in = Fout * (Cc + Ce);
%     
%     Hcin = Hc;
%     Hein = He;
%     Hcout = CpC * (Tcond - Ta);
%     Heout = CpE * (Tcond - Ta);
%     Qcond = UAcond*(Tcond - Ta);
%     
%     lamdac = 0.3445;
%     lamdae = 0.1256;
%     
%     
%     dyc_dt = 0;
%     dye_dt = 0;
%     
%     if NC > 0 && NE > 0
%         den_n_liq = (yc*lamdac + ye*lamdae) - (yc*Hcout + ye*Heout);
%         n_liq = (n_C_in*Hcin + n_E_in*Hein - (n_CE_in*(yc*Hcout+ye*Heout)) - Qcond)/den_n_liq;
%         dyc_dt = (n_C_in - yc*(n_CE_in - (2 * n_liq)))/NCE;
%         dye_dt = (n_E_in - ye*(n_CE_in - (2 * n_liq)))/NCE;
%     else
%         dyc_dt = 0;
%         dye_dt = 0;
%     end
%     
%     dyc_dt = n_C_in/NCE;
%     dye_dt = n_E_in/NCE;
    
    
    Y = [dCa_dt; dCb_dt; dCc_dt; dCd_dt; dCe_dt; dT1_dt];
end