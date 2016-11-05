function [a, omega] = true_acc_vel(t)
    %% Acceleration
    a_imu = sin(t);
    w_1 = 0.4 * wgn(1,length(t),1);
    
    T_1 = 100;
    w_2 = 0.05 * wgn(1,length(t),1);
    b_1 = zeros(1,length(t));
    b_1(1) = w_2(1);
    b_1_dot = zeros(length(t));
    
    for i = 1:length(t)-1
        b_1_dot(i) = -(1/T_1)*b_1(i) + w_2(i);
        b_1(i+1) = b_1(i) + b_1_dot(i);
    end
    
    a = a_imu + w_1 + b_1;
    
    %% Angular velocity
    omega_imu = sin(2*t);
    w_3 = 0.4 * wgn(1,length(t),1);
    
    T_2 = 100;
    w_4 = 0.05 * wgn(1,length(t),1);
    b_2 = zeros(1,length(t));
    b_2(1) = w_4(1);
    b_2_dot = zeros(length(t));
    
    for i = 1:length(t)-1
        b_2_dot(i) = -(1/T_2)*b_2(i) + w_4(i);
        b_2(i+1) = b_2(i) + b_2_dot(i);
    end
    
    omega = omega_imu + w_3 + b_2;

end