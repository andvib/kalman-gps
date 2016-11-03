function [a, omega] = true_acc_vel(t)
    %% Acceleration
    a_imu = sin(t);
    w_1 = 0.5 * wgn(1,length(t),1);
    
    a = a_imu + w_1;
    
    %% Angular velocity
    omega_imu = sin(t);
    w_2 = 0.5 * wgn(1,length(t),1);
    
    omega = omega_imu + w_2;

end