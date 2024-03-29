clear all; close all; clc;

%% variables
time_step = 0.001;                                                          %time steps (s)
time_end = 1;                                                               %length of model (s)
U = 4.5;                                                                    %Voltage (V)   

%% Constants
m_total = 0.7179;                                                           %total vibrating mass (kg)
m_bolt_prop = 0.0163;                                                       %mass of bolt and the propeller (kg)
l = 0.0061;                                                                 %lenth between center of rotation and center of mass (m)
f_n = 15;                                                                   %natural frequency (Hz)
omega_n = f_n * 2 * pi;                                                     %natural frequency (rad/s)
k = m_total * omega_n.^2;                                                   %spring constant (N/m)
zeta = 0.019;                                                               %damping factor (-)
g = 9.81;                                                                   %gravitational acceleration (m/s^2)
F_g = g*m_bolt_prop;                                                        %gravitational force (N)

%% Determining Angle(theta) and Angular velocity(omega) over angular acceleration(alpha), based on the gravitational force(F_g)
n = 1;
theta(n,1) = 0;
theta(n,2) = 0;
omega(n,1) = 0;
omega_av = -0.1706*U.^2 + 19.5596*U - 18.8967;                              %experimentaly determined average angular velocity (rad/s) 
omega(n,2) = omega_av;                                                      %for now, this maybe have to be calculated later
for i = 0:time_step:time_end;
    n = n+1;
    theta(n,1) = i;
    omega(n,1) = i;
    alpha(n,1) = i;
    alpha(n,2) = (-cos(theta(n-1,2)) * g)./(l);                             %angular acceleration of the propeller (rad/s^2)
    omega(n,2) = omega(n-1,2) + alpha(n,2)*time_step;                       %angular velocity of the propeller (rad/s)
    theta(n,2) = theta(n-1,2) + omega(n,2)*time_step;                       %angle of the propellor (rad)
end

%% Centripetal force(F_c) over Angular velocity(omega)
F_c = (m_bolt_prop *omega_av.^2*l);

%% Amplitude(X) over Angular velocity(omega), Centripetal force(F_c)
n = 1;
for i = 0:time_step:time_end
    n = n+1;
    X(n,1) = i;
    X(n,2) = (F_c./k)./(((1-(omega_av/omega_n).^2).^2 + (2 * zeta * (omega_av./omega_n)).^2)^0.5);       %amplitude(m), dynamics book
    X(n,3) = X(n,2).*cos(theta(n,2));                                                                    %x(t)  = displacement over time 
end
%% Frequency(f) over Angular velocity(omega)
n = 1;
for i = 0:time_step:time_end
    n = n+1;
    f(n,1) = i;
    f(n,2) = (omega(n,2))./(2*pi);                                          %frequency (Hz)
end

%% for loop centripetal force 0->4V
  n=0;
  p=0;
for j= 0:0.01:0.98
   n=n+1;
   Range_U(n,1) = j;
   Range_omega(n,1) = 0;      %omega, angular velocity
   Fc_omega(n,1)= m_bolt_prop*Range_omega(n,1).^2*l;                                 %centripedal force
end
%% for loop centripetal force 4->24V
for j= 0.99:0.01:24
   n=n+1;
   Range_U(n,1) = j;
   Range_omega(n,1) = -0.1706*Range_U(n,1).^2 + 19.5596*Range_U(n,1) - 18.8967;      %omega, angular velocity
   Fc_omega(n,1)= m_bolt_prop*Range_omega(n,1).^2*l;                                 %centripedal force
end
%% Resonance centipedal force

X_omega(:,1) = (Fc_omega(:,1)./k)./(((1 - (Range_omega(:,1)./omega_n).^2).^2 + (2 * zeta * (Range_omega(:,1)./omega_n)).^2).^0.5);       %formula for excitation from dynamics book













%% plots
 figure(1);
 yyaxis left
 plot (X(:,1),X(:,3))
 xlabel('Time (s)')
 ylabel('Amplitude (m)')
 title('Distance from equilibrium over time')
 
 yyaxis right
 plot(omega(:,1), omega(:,2))
 xlabel('Time (s)')
 ylabel('Angular velocity \omega of the propeller (rad/s)')
 title('Angular velocity over time')
 
 figure(2);
 plot(theta(:,1), theta(:,2))
 xlabel('Time (s)')
 ylabel('Angle theta of the propeller (rad)')
 title('Angle of the propeller over time')

figure(3);
steps = 349/2300;
angular_velocity= [0:steps:349];
hAX=axes;                 % first axes, save handle
pos=get(hAX,'position')   % get the position vector
pos = [ 0.1300    0.1100    0.7750    0.8150];
pos1=pos(2);              % save the original bottom position
pos(2)=pos(2)+pos1; pos(4)=pos(4)-pos1;  % raise bottom/reduce height->same overall upper position
set(hAX,'position',pos)   % and resize first axes
pos(2)=pos1; pos(4)=0.01; % reset bottom to original and small height
plot(Range_U(:,1), X_omega(:,1))
title('Displacement of the top of the system')

hAX(2)=axes('position',pos,'color','none');  % and create the second
plot(angular_velocity, 0)
ylabel(hAX(1),'Displacement of system [m]')
xlabel(hAX(1),'Supplied voltage [V]')
xlabel(hAX(2),'Angular velocity [rad/s]')
set(hAX(2),'xcolor','r','ycolor','r')
