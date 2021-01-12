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

end