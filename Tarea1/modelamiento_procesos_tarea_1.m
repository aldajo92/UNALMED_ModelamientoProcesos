clc; clear; close all;
% Parameters
p.Ro = 1.225;       %[kg/m3]
p.V = 160;          %[m3]
p.Kv1 = 0.315;      %[m2]
p.Kv2 = 0.25;       %[m2]
p.Kv3 = 0.35;       %[m2]
p.P1 = 25;          %[bar]
p.P2 = 20;          %[bar]
p.Patm = 0.8535;    %[bar][Medellin]
% initial condition
x0 = [4];

tv = [0 600];
[t, result] = ode45(@dynamic, tv, x0, [], p);
plot(t, result);

function dP = dynamic(t, x, p)
    P = x;
    e1 = p.Kv1*sqrt((p.P1 - P)/p.Ro);
    e2 = p.Kv2*sqrt((p.P2 - P)/p.Ro);
    e3 = p.Kv3*sqrt((P - p.Patm)/p.Ro);
    dP = P*(e1 + e2 - e3)/p.V;
end

% % parameters
% p.At = 19.7;        %[m2]
% p.Kv1= 0.769;       %[m2]
% p.Kv2 = 0.694;      %[m2]
% p.Fi = 7.2;         %[m3/s]
% p.g = 9.8;          %[m/s2]
% % initial condition
% h1 = 4;
% h2 = 2;
% x0 = [h1; h2];
% % time
% tv = [0 600];
% [t, h] = ode45(@dinamica, tv, x0, [], p);
% plot(t, h);
% 
% function dh = dinamica(t, x, p)
%     h1 = x(1);
%     h2 = x(2);
%     dh1dt = (p.Fi-p.Kv1*sqrt(p.g*h1))/p.At;
%     dh2dt = (p.Kv1*sqrt(p.g*h1)-p.Kv2*sqrt(p.g*h2))/p.At;
%     dh = [dh1dt; dh2dt];
% end