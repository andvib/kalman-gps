close all;
clear all;

load('data/data.mat');

% Set initial position estimate
x0 = [0, 0, 0, 0]';
x = x0;
x_log = [];


% Loop through all measurements
for t = 1:7200
    visible_pos = [];
    visible_pseudo = [];
    visible_el = [];
    P_delta_log = [];
    
    % Loop through measurements to find visible satellites
    for j = 1:32
        if Satpos(:,j,t) ~= 0;
            visible_pos = [visible_pos, Satpos(:,j,t)];
            visible_pseudo = [visible_pseudo, PR(j, t)];
            visible_el = [visible_el, EL(j,t)];
        end
    end
    no_visible = length(visible_pos);
    
    % Estimation loop
    for i = 1:10
        % Estimate pseudoranges
        P_hat = sqrt((x(1) - visible_pos(1,:)).^2 + ...
                     (x(2) - visible_pos(2,:)).^2 + ...
                     (x(3) - visible_pos(3,:)).^2);

        % Observed minus computed ranges
        P_delta = P_hat - visible_pseudo;
        P_delta_log = [P_delta_log, P_delta'];
        
        % Geometry matrix
        G = [(visible_pos(1,:)' - x(1))./P_hat(:), ...
             (visible_pos(2,:)' - x(2))./P_hat(:), ...
             (visible_pos(3,:)' - x(3))./P_hat(:), ...
             ones(no_visible,1)];

        % Calculate covariance matrix
        %std_deviation = std(P_delta_log(:,:));
        std_deviation = 20;
        R = diag(std_deviation.^2 .* (1./(sind(visible_el).^2)));
         
        % Determine the offset delta_x
        %x_delta = inv(G'*inv(R)*G) * G' * inv(R) * P_delta';
        x_delta = inv(G'*G) * G' * P_delta';
        x = x + x_delta;
    end
    x_log = [x_log, x(1:3)];
end

% Transform from ECEF to ENU
wgs84 = wgs84Ellipsoid;

lla_pos = ecef2lla(P0');
[xP0, yP0, zP0] = ecef2enu(P0(1), P0(2), P0(3), ...
                           lla_pos(1), lla_pos(2), lla_pos(3), wgs84);

lla = ecef2lla(x_log')';
[xEast, yNorth, zUp] = ecef2enu(x_log(1,:),x_log(2,:),x_log(3,:), ...
                                lla_pos(1), lla_pos(2), lla_pos(3), wgs84);

                       
% Plot position error in East, North and Vertical frame
figure(1);
plot((xEast - xP0));
grid on;
title('Position error, east')

figure(2);
plot((yNorth - yP0));
grid on;
title('Position error, north')

figure(3);
plot((zUp - zP0));
grid on;
title('Position error, up')

figure(4)
scatter(xP0,yP0,'b')
hold on
scatter(xEast(7200),yNorth(7200),'r')