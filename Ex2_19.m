%% 
% Ex 2.19 Problem statement 

% Plot the variations of the natural frequency and the time period with 
% static deflection of an undamped system using MATLAB.

% Code modified from Mechanical Vibrations - S.S. Rao
clc
close all
clear
g = 9.81;                                   %m/s^2
for i = 1: 101                              
t(i) = 0.01 + (0.5-0.01) * (i-1)/100;
w(i) = (g/t(i))^0.5;
tao(i) = 2 * pi * (t(i)/g)^0.5;
end
plot(t,w);
txt = {'\omega_n'}
text(0.02,25,txt,'FontSize',14,'Color','black')
hold on;
plot(t, tao);
txt = {'T_n'}
text(0.05,2,txt,'FontSize',12,'Color','black')
xlabel('\delta_s_t');
ylabel('
title('Example 2.19');
