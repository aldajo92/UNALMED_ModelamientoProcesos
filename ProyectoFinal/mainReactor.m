clc; clear; close all;
%% Tiempo de simulacion
time = [0:0.001:10];


%% Condiciones iniciales
% X0 = [CA;CB;CC;CD;CE;CF;T;Fin];
X0 = [1;2;3;4;5;6;7;8;9];

%% Simulacion
reactor(time, X0)
% [t, result] = ode45(@horno, time, X0, [], params);

