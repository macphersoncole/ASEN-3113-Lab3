close all;
clear;
clc;

T_max = 30 + 273.15;
T_min = 20 + 273.15;
T_standby = -40 + 273.15;

Q_sat = 20;
e = 0.85;
a = 0.2;
sigma = 5.670373e-8;
G_s = 1361;
A = 0.286448;
q_eq = 75.5;
q_ec = 11;

time = 0:(24*60); % minutes in a day
theta = linspace(0,360,length(time)) - 90;
Q_heater = zeros(1,length(time));
Q_environment = zeros(1,length(time));
Q_sat_vec = zeros(1,length(time)) + 20;
T = zeros(1,length(time));

%% In operation

for i = 1:length(time)
    
    if theta(i) < (90-atan(1/6.6)) || theta(i) >= 270 % Sun hitting radiator
        T_temp = ((Q_sat+e*q_eq*A+G_s*A*a*cosd(theta(i)))/(e*sigma*A))^(1/4);
        Q_environment(i) = e*q_eq*A+G_s*A*a*cosd(theta(i));
        if T_temp < T_min
            Q_heater(i) = e*sigma*A*T_min^(4)-Q_sat-e*q_eq*A-G_s*A*a*cosd(theta(i));
            T(i) = T_min;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    elseif theta(i) <= 90 && theta(i) >= (90-atan(1/6.6)) % Eclipsed reaitaor
        T_temp = ((Q_sat+e*q_ec*A)/(e*sigma*A))^(1/4);
        Q_environment(i) = e*q_eq*A;
        if T_temp < T_min
            Q_heater(i) = e*sigma*A*T_min^(4)-e*q_ec*A-Q_sat;
            T(i) = T_min;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    elseif theta(i) < 270 && theta(i) > 90 % No sun hitting the radiator
        T_temp = ((Q_sat)/(e*sigma*A))^(1/4);
        if T_temp < T_min
            Q_heater(i) = e*sigma*A*T_min^(4)-Q_sat;
            T(i) = T_min;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    end
    
end

figure;
subplot(2,1,1);
plot(time,Q_heater,'r',time,Q_sat_vec,'b',time,Q_environment,'g');
grid on;
title('Q_{in} vs Time - Equinox - Instrument on - Survival Temp');
xlabel('Time [minutes]');
ylabel('Heater Input [W]');
xlim([0 time(end)]);
legend('Q_{in} Spacecraft Heater','Q_{in} Instrument','Q_{in} Environment','location','east');

subplot(2,1,2);
plot(time,T-273.15,'r');
grid on;
title('Temperature vs Time - Equinox - Instrument on - Survival Temp');
xlabel('Time [minutes]');
ylabel('Temperature [^{o}C]');
xlim([0 time(end)]);

%% In Standby

for i = 1:length(time)
    
    if theta(i) < (90-atan(1/6.6)) || theta(i) >= 270 % Sun hitting radiator
        T_temp = ((e*q_eq*A+G_s*A*a*cosd(theta(i)))/(e*sigma*A))^(1/4);
        Q_environment(i) = e*q_eq*A+G_s*A*a*cosd(theta(i));
        if T_temp < T_standby
            Q_heater(i) = e*sigma*A*T_standby^(4)-e*q_eq*A-G_s*A*a*cosd(theta(i));
            T(i) = T_standby;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    elseif theta(i) <= 90 && theta(i) >= (90-atan(1/6.6)) % Eclipsed reaitaor
        T_temp = ((e*q_ec*A)/(e*sigma*A))^(1/4);
        Q_environment(i) = e*q_eq*A;
        if T_temp < T_standby
            Q_heater(i) = e*sigma*A*T_standby^(4)-e*q_ec*A;
            T(i) = T_standby;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    elseif theta(i) < 270 && theta(i) > 90 % No sun hitting the radiator
        T_temp = ((0)/(e*sigma*A))^(1/4);
        if T_temp < T_standby
            Q_heater(i) = e*sigma*A*T_standby^(4);
            T(i) = T_standby;
        else
            Q_heater(i) = 0;
            T(i) = T_temp;
        end
    end
    
end

figure;
subplot(2,1,1);
plot(time,Q_heater,'r',time,Q_sat_vec,'b',time,Q_environment,'g');
grid on;
title('Q_{in} vs Time - Equinox - Instrument on - Survival Temp');
xlabel('Time [minutes]');
ylabel('Heater Input [W]');
xlim([0 time(end)]);
legend('Q_{in} Spacecraft Heater','Q_{in} Instrument','Q_{in} Environment','location','northeast');

subplot(2,1,2);
plot(time,T-273.15,'r');
grid on;
title('Temperature vs Time - Equinox - Instrument on - Survival Temp');
xlabel('Time [minutes]');
ylabel('Temperature [^{o}C]');
xlim([0 time(end)]);

