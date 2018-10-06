
%% MEM 423 - Mechanics of Vibration - Fall 2018
% Example 2.22) Develop a general-purpose MATLAB program, to find the free-vibration
%               response of a viscously damped system. Use the program to find the response of a system with the
%               following data:

clc
clear all
close

%% Data:
m= 450;     % mass
k= 26519.2; % spring stiffness
c= 1000;    % damping constant   
x0= 5.3965700e-01;    % initial displacement
xd0= 1.0;   % initial velocity

% Discretization
n= 100;     %total no. of increments
delta= 2.50e-2; %increment size
t = [0: delta: delta*n];

%Data Initilization
x=zeros(1,101);
xd=zeros(1,101);
xdd=zeros(1,101);
x(1)=x0;
xd(1)=xd0;

w_n=sqrt(k/m);  %Undamped natural frequency
zeta = c/(2*m*w_n); %Damping co-efficient


if zeta<0
    disp ('check input parameters')

elseif zeta==0
    disp('system is undamped')
    
elseif zeta==1
    disp ('system is critically damped')
    for i=1:(2.5/delta)
        x(i)=(x0+(xd0+w_n*x0)*t(i))*exp(-w_n*t(i));
    end
    
elseif zeta>1
    disp ('system is over damped')
    c1=(x0*w_n*(zeta+sqrt((zeta^2)-1))+xd0)/(2*w_n*sqrt((zeta^2)-1));
    c2=(-x0*w_n*(-zeta-sqrt((zeta^2)-1))-xd0)/(2*w_n*sqrt((zeta^2)-1));
    for i=1:(2.5/delta)
        x(i+1)=(c1*exp((-zeta+sqrt((zeta^2-1)))*w_n*t(i)))+(c2*exp((-zeta-sqrt((zeta^2-1)))*w_n*t(i)));
    end
    
else
    disp ('system is under damped')
    w_d=(sqrt(1-(zeta^2)))*w_n      %Damped natural frequency
    bound(1)=(sqrt(x0^2*w_n^2+xd0^2+2*x0+xd0*zeta*w_n)/w_d);    
    for i=2:(2.5/delta)+1
        x(i)=(exp(-zeta*w_n*t(i)))*((x0*cos(w_d*t(i)))+(((xd0+(zeta*w_n*x0))/w_d)*sin(w_d*t(i))));  % from week 2 notes - Eq (7) | Eq 2.72 in textbook
        xd(i)=(-zeta*w_n*x(i))+(exp(-zeta*w_n*t(i)))*((-x0*w_d*sin(w_d*t(i)))+((xd0+(zeta*w_n*x0))*cos(w_d*t(i))));
        xdd(i)=(-zeta*w_n*xd(i))-(zeta*w_n*exp(-zeta*w_n*t(i)))*((-x0*w_d*sin(w_d*t(i)))+((xd0+(zeta*w_n*x0))*cos(w_d*t(i))))-(w_d^2*x(i));
        bound(i)=(sqrt(x0^2*w_n^2+xd0^2+2*x0+xd0*zeta*w_n)/w_d)*exp(-zeta*w_n*t(i));
    end
    xdd(1)=(-zeta*w_n*xd(1))-(zeta*w_n*exp(-zeta*w_n*t(1)))*((-x0*w_d*sin(w_d*t(1)))+((xd0+(zeta*w_n*x0))*cos(w_d*t(i))))-(w_d^2*x(1));
end

%% Plot
figure(1)
plot(t,x,'r');
hold on
plot(t,bound,'k--');
plot(t,-bound,'k--');
hold off
xlabel('Time, t');
ylabel('Displacement, x');
title('Under damped vibration - displacement');
legend('Displacement','upper bound','lower bound');

figure(2)
plot(t,x);
hold on
plot(t,xd,'--');
plot(t,xdd,'o-');
xlabel('Time, t');
ylabel('Displacement, Velocity, Acceleration');
title('Under damped vibration');
legend('Displacement','Velocity','Acceleration');