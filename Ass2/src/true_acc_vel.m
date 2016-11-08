function [a, omega] = true_acc_vel(t)
    %% Acceleration
    a = 0.5*sin(t);
    
    %% Angular velocity
    omega = 0.5*sin(t);

end