%ASEN 3113 Thermo Lab 3 Design Lab
% Nate Kuczun 

%Process
%{
- find max Qin to radiator at any time
- size radiator/ material 
- find Qin at a specific time based on sunlight, instrument, spacecraft IR
- calculate heater to compensate for Qin loss
- plot results for solstices and equinox 

%}


%% - Setup
    clear 
    close all
    clc

    %dimesions (starting point) 
    radLength = 1; %m
    radWidth = 1; %m
    radArea = radLength*radWidth; %m^2

    %temp ratings
    tempMax_op = 30 + 273; %Kelvin
    tempMin_op = 20 + 273; %Kelvin
    tempMin_stb = -40 + 273; %Standby temp (in kelvin) 

    sigma = 5.670 * 10^(-8); % W/m2-K4
    Gs = 1361; %w/m^2 


    %surface
    alpha = 0.2; 
    %{
    %highly pos alum 
    alpha = 0.1 - 0.4
    %rough alum
    alpha = 0.4-0.65

    %}
    eps = 0.85; %emissivity 
    %{
    %highly pos alum 
    e = 0.039-0.057
    %ann alum
    e = 0.77; 
    %rough alum
    e = 0.07; 
    %}

%% Find Raditor pannel size
    thetaWinter = 23.5; %degrees
    thetaForMaxSun = 0; %degrees

    fprintf('Sizing radiator ...\n')
    error = 1; 
    j = 1; 
    while error > 0.01
        %find max heat generated
        Qin_inst = 20; %20W continously dumped into radiator from instrument
        qin_spacecraft_max = (88 + 63)/2; %w/m^2 during equ.  
        Qin_spacecraft_max = eps * qin_spacecraft_max * radArea; %w  
        Qin_sunlight_max = Gs * radArea * alpha * cosd(thetaForMaxSun); %w  

        Qin_env_max = Qin_inst + Qin_spacecraft_max + Qin_sunlight_max; 


        %% Size radiator/ material 
        Qout = Qin_env_max;
        %assume radiator will always be in middle of thermal range
        T = tempMax_op;  
        radAreaNew = Qout / ( eps * sigma * T^4 ); %m^2 
        error = abs(radArea - radAreaNew)/radAreaNew * 100; 
        radArea =  radArea + (radAreaNew-radArea)/2; 

        j = j + 1; 
    end

    fprintf('%d sizing iterations. Final area %d m^2 \n',j, radArea)

%% Time vector 

    [~,~,time] = radiatorAngles(); %get radiator angle time base
    t = time; 

%% Winter solstice 

% spacecraft instrument is on 
% IR backload = 88 w/m^2
% spacecraft orbit (not radiator) is in continual sunlight 

    active = 'y'; 
    t = 0:(24*60); %one day 
    
    %Heat transfer from enviroment 
    for i = 1:length(t)
        Qin_env(i) = QinFunctionWinter(t(i),radArea,Gs,alpha,eps);
    end

    %Calculate heater to compensate for Qin loss
    Qin_heater = Qin_env_max - Qin_env - Qin_inst; %Electronic heater - compensate for max heat

    %find temperature 
    Qin_total = Qin_heater + Qin_env + Qin_inst; 
    Qout = Qin_total;
    T_spacecraft = (Qout ./ (radArea .* sigma .* eps)).^(1/4); %K 
    T_spacecraft = T_spacecraft - 273; %convert to C 
    
    %plot 
    figure(1)
    
    subplot(2,1,1)
    plot(t,Qin_env,'g-')
    hold on;
    yline(Qin_inst,'b-');
    plot(t,Qin_heater,'r-')
    xlabel('Time [Minutes]')
    ylabel('Q_{in} [W]');
    title('Q_{in} vs Time - Winter Solstice - Instrument On - Operational Temp');
    legend('Q_{in} Enviroment','Q_{in} Instrument','Q_{in} Spacecraft heater'); 
    set(legend, 'Location', 'Best')
    hold off; 
    
    subplot(2,1,2)
    plot(t,T_spacecraft,'r-')
    xlabel('Time [Minutes]')
    ylabel('Temperature [C]');
    title('Q_{in} vs Time - Winter Solstice - Instrument off - Survival Temp');
    
% spacecraft instrument is off 
% IR backload = 88 w/m^2
% spacecraft orbit (not radiator) is in continual sunlight 

    active = 'n'; 
    t = 0:(24*60); %one day 
    
    T = tempMin_stb; %survival temp 
    Qin_env_max_survial = radArea * eps *sigma * T^4; %find Qout radiated at survival temp 
    
    %Heat transfer from enviroment 
    for i = 1:length(t)
        Qin_env(i) = QinFunctionWinter(t(i),radArea,Gs,alpha,eps); 
    end

    %Calculate heater to compensate for Qin loss
    Qin_heater = Qin_env_max_survial - Qin_env; %Electronic heater - compensate for max heat 
    Qin_heater = max(Qin_heater,0); %set negative heater values to zero 
    
    %find temperature 
    Qin_total = Qin_heater + Qin_env; 
    Qout = Qin_total;
    T_spacecraft = (Qout ./ (radArea .* sigma .* eps)).^(1/4); %K 
    T_spacecraft = T_spacecraft - 273; %convert to C 

    %plot 
    figure(2)
    
    subplot(2,1,1)
    plot(t,Qin_env,'g-')
    hold on;
    yline(Qin_inst,'b-');
    plot(t,Qin_heater,'r-')
    xlabel('Time [Minutes]')
    ylabel('Q_{in} [W]');
    title('Q_{in} vs Time - Winter Solstice - Instrument off - Survival Temp');
    legend('Q_{in} Enviroment','Q_{in} Instrument','Q_{in} Spacecraft heater'); 
    set(legend, 'Location', 'Best')
    hold off; 
    
    subplot(2,1,2)
    plot(t,T_spacecraft,'r-')
    xlabel('Time [Minutes]')
    ylabel('Temperature [C]');
    title('Q_{in} vs Time - Winter Solstice - Instrument off - Survival Temp');

    
%% Summer solstice 

% spacecraft instrument is on 
% IR backload = 63 w/m^2
% spacecraft orbit (not radiator) is in continual sunlight 

active = 'y'; 
    t = 0:(24*60); %one day 
    
    %Heat transfer from enviroment 
    for i = 1:length(t)
        Qin_env(i) = QinFunctionSummer(t(i),radArea,Gs,alpha,eps);
    end

    %Calculate heater to compensate for Qin loss
    Qin_heater = Qin_env_max - Qin_env - Qin_inst; %Electronic heater - compensate for max heat

    %find temperature 
    Qin_total = Qin_heater + Qin_env + Qin_inst; 
    Qout = Qin_total;
    T_spacecraft = (Qout ./ (radArea .* sigma .* eps)).^(1/4); %K 
    T_spacecraft = T_spacecraft - 273; %convert to C 
    
    %plot 
    figure(3)
    
    subplot(2,1,1)
    plot(t,Qin_env,'g-')
    hold on;
    yline(Qin_inst,'b-');
    plot(t,Qin_heater,'r-')
    xlabel('Time [Minutes]')
    ylabel('Q_{in} [W]');
    title('Q_{in} vs Time - Summer Solstice - Instrument On - Operational Temp');
    legend('Q_{in} Enviroment','Q_{in} Instrument','Q_{in} Spacecraft heater');
    set(legend, 'Location', 'Best')
    hold off; 
    
    subplot(2,1,2)
    plot(t,T_spacecraft,'r-')
    xlabel('Time [Minutes]')
    ylabel('Temperature [C]');
    title('Q_{in} vs Time - Summer Solstice - Instrument off - Survival Temp');
    
% spacecraft instrument is off 
% IR backload = 88 w/m^2
% spacecraft orbit (not radiator) is in continual sunlight 

    active = 'n'; 
    t = 0:(24*60); %one day 
    
    T = tempMin_stb; %survival temp 
    Qin_env_max_survial = radArea * eps *sigma * T^4; %find Qout radiated at survival temp 
    
    %Heat transfer from enviroment 
    for i = 1:length(t)
        Qin_env(i) = QinFunctionSummer(t(i),radArea,Gs,alpha,eps); 
    end

    %Calculate heater to compensate for Qin loss
    Qin_heater = Qin_env_max_survial - Qin_env; %Electronic heater - compensate for max heat 
    Qin_heater = max(Qin_heater,0); %set negative heater values to zero 
    
    %find temperature 
    Qin_total = Qin_heater + Qin_env; 
    Qout = Qin_total;
    T_spacecraft = (Qout ./ (radArea .* sigma .* eps)).^(1/4); %K 
    T_spacecraft = T_spacecraft - 273; %convert to C 

    %plot 
    figure(4)
    
    subplot(2,1,1)
    plot(t,Qin_env,'g-')
    hold on;
    yline(Qin_inst,'b-');
    plot(t,Qin_heater,'r-')
    xlabel('Time [Minutes]')
    ylabel('Q_{in} [W]');
    title('Q_{in} vs Time - Summer Solstice - Instrument off - Survival Temp');
    legend('Q_{in} Enviroment','Q_{in} Instrument','Q_{in} Spacecraft heater');
    set(legend, 'Location', 'Best')
    hold off; 
    
    subplot(2,1,2)
    plot(t,T_spacecraft,'r-')
    xlabel('Time [Minutes]')
    ylabel('Temperature [C]');
    title('Q_{in} vs Time - Summer Solstice - Instrument off - Survival Temp');

    
%% FUNCTIONS

%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%

function[Qin_env] = QinFunctionWinter(t,radArea,Gs,alpha,eps)

qin_spacecraft = 88; %w/m^2
Qin_spacecraft = eps * radArea*qin_spacecraft; %w

Qin_sunlight = Gs * radArea * alpha * cosd(23.5)*cosd(thetaFuncWinter(t)); %w
if Qin_sunlight < 0 
    Qin_sunlight = 0;
end 
Qin_env = Qin_spacecraft + Qin_sunlight; 

end

function[Qin_env] = QinFunctionSummer(t,radArea,Gs,alpha,eps)

qin_spacecraft = 63; %w/m^2
Qin_spacecraft = eps * radArea*qin_spacecraft; %w

Qin_sunlight = Gs * radArea * alpha * cosd(thetaFuncSummer(t)); %w
if Qin_sunlight < 0 
    Qin_sunlight = 0;
end 
Qin_env = Qin_spacecraft + Qin_sunlight; 

end

function[theta] = thetaFuncWinter(t)
%input t = time [min] 
%output theta [degrees] 

%spacecraft is always in sunlight
%spacecraft is pointed at earth
%radiator is in the back of the spacecraft 
%in northern hemisphere winter earth tilts away from sun --> radiator gets
%sunlight at noon 
theta = -90 + 90*sin(2*pi*t/(60*24)); %perpendicular off set + angle relative to sun starting at noon

end


function[theta] = thetaFuncSummer(t)
%input t = time [min] 
%output theta [degrees] 

%spacecraft is always in sunlight
%spacecraft is pointed at earth
%radiator is in the back of the spacecraft 
%in northern hemisphere summer earth tilts toward sun --> radiator gets
%sunlight at midnight 
theta = -90 + 90*sin(-2*pi*t/(60*24)); %earth tilt + perpendicular off set + angle relative to sun starting at noon

end

%{
function[Qin_env] = qinFunction(t,radArea,Gs,alpha)

%assume starts in the winter (GOES R launched Nov 2016) 
%assume that spacecraft acts as thermal res. 
qin_spacecraft = @(t) 12.5 * cos( (2*pi) / (365*24) *t) + 75.5; %function of time over a year W/m^2
Qin_craft_day = radArea * qin_spacecraft(t); %W 
Qin_craft_night = radArea * 11; %w/m^2 

%theta function over time
thetaFunc = @(t) t; %degrees 

%sunlight
Qin_sunlight = Gs * radArea * alpha * cosd(thetaFunc(t)); %w 

%eclipse function
eclipseFunc = @(t) round(sin(2.* pi ./ (24 ).* t)); % change later

for i = t(1):t(end)
    if eclipseFunc(t) == 1
        Qin_env = Qin_craft_day + Qin_sunlight;
    else
        Qin_env = Qin_craft_night;
    end 
end

end
%}

function [theta_true,theta_unchanged,time] = radiatorAngles()

    time = 0:(1/24):365.25; % days
    theta_unchanged = linspace(0,360*365.25,length(time)); % theta of the radiator as it rotates around earth
    j = 1;
    for i = 1:24:length(time)
        if i == 1
        elseif i == 8761
            theta_unchanged(i:(i+6)) = theta_unchanged(i:(i+6)) - (360*j);
        else
            theta_unchanged(i:(i+23)) = theta_unchanged(i:(i+23)) - (360*j);
            j = j+1;
        end
    end
    
    theta_true = zeros(1,length(time));
    for i = 1:length(theta_unchanged)
        if cosd(theta_unchanged(i)) < 0
            theta_true(i) = 450;
        elseif cosd(theta_unchanged(i)) >= 0
            theta_true(i) = theta_unchanged(i);
        end
    end
    
    %{
    figure
    plot(time,cosd(theta_unchanged));
    title('Unchanged Angles (w/ negative values)');
    xlabel('Time [daya]');
    ylabel('cosd(thata)');
    
    figure
    plot(time,cosd(theta_true));
    title('Changed Angles (accounts for time where there is no sun)');
    xlabel('Time [daya]');
    ylabel('cosd(thata)');
    %}
    
end