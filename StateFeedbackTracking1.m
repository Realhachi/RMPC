%state feedback controller to stablize output
clc;clear;close all;
% load Control Package
% pkg load control
%% define sim span
tspan = [0,40];
z0 = [0 0];
%% Solve 
[t, z] = ode45(@sys,tspan,z0);
hold on;
plot (t,z(:,2));
% hold on
% plot(t,u)
grid on;