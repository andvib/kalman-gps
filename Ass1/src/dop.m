close all;
clear all;

load('data/data.mat');

GDOP = [];
PDOP = [];
HDOP = [];
VDOP = [];
TDOP = [];
no_satellits = [];

for t = 1:length(PR)
    visible_pos = [];
    
    % Loop through measurements to find visible satellites
    for j = 1:32
        if Satpos(:,j,t) ~= 0;
            visible_pos = [visible_pos, Satpos(:,j,t)];
        end
    end
    no_satellits = [no_satellits length(visible_pos)];
    % Calculate unit vectors
    R = sqrt((P0(1) - visible_pos(1,:)).^2 + ...
             (P0(2) - visible_pos(2,:)).^2 + ...
             (P0(3) - visible_pos(3,:)).^2);
    h_x = (visible_pos(1,:) - P0(1))./R(:)';
    h_y = (visible_pos(2,:) - P0(2))./R(:)';
    h_z = (visible_pos(3,:) - P0(3))./R(:)';
    
    G = [h_x' h_y' h_z' -ones(length(visible_pos),1)];
    
    H = inv(G' * G);
    
    % Calculate dilution of precision
    gdop = sqrt(H(1,1) + H(2,2) + H(3,3) + H(4,4));
    pdop = sqrt(H(1,1) + H(2,2) + H(3,3));
    hdop = sqrt(H(1,1) + H(2,2));
    vdop = sqrt(H(3,3));
    tdop = sqrt(H(4,4));
    
    % Store the results
    GDOP = [GDOP gdop];
    PDOP = [PDOP pdop];
    HDOP = [HDOP hdop];
    VDOP = [VDOP vdop];
    TDOP = [TDOP tdop];
end

% Plotting
figure(1);
hold on;
plot(GDOP);
plot(PDOP);
plot(HDOP);
plot(VDOP);
plot(TDOP);
legend('GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP', ...
       'Location','northeastoutside');
grid on;
saveas(gcf, 'dop', 'epsc')

figure(2)
plot(no_satellits);
grid on;
saveas(gcf, 'dop_no_satellites', 'epsc');