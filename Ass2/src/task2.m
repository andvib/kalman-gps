close all;
clear all;

t = [0:0.01:50];

global T_1 T_2;
T_1 = 1;
T_2 = 1;

%% Generate acceleration and angular velocity
[a, omega] = true_acc_vel(t);

figure(1);
subplot(2,1,1);
plot(t, a);
title('True acceleration a');
grid on;

subplot(2,1,2);
plot(t, omega);
title('Angular velocity \omega');
grid on;

saveas(gcf, 'task1', 'epsc');

%% Simulate discretized system
u = [a; omega];

x0 = [0;0;0];

x = zeros(3,length(t));
x(:,1) = x0;

states = disc_sys(x, u, t);

figure(2);
subplot(3,1,1)
plot(t, states(1,:));
grid on;
ylabel('[m]');
title('Position x')

subplot(3,1,2);
plot(t, states(2,:));
grid on;
ylabel('[m/s]');
title('Velocity v')

subplot(3,1,3);
plot(t, states(3,:));
grid on;
ylabel('[rad/s]'),
xlabel('[s]')
title('Angular velocity \theta')

saveas(gcf, 'disc_states', 'epsc');

%% Simulate discrete bias and white noise
[bias, white, vhite] = disc_bias(t);

figure(3);
plot(t, bias);
legend('b_1', 'b_2');
xlabel('[t]');
grid on;
saveas(gcf, 'disc_bias', 'epsc');

%% Simulate system with noise and bias
a_sig = a + bias(1,:) + white(1,:);
omega_sig = omega + bias(2,:) + white(2,:);

states_sig = zeros(3,length(t));
states_sig = disc_sys(states_sig, [a_sig ; omega_sig], t);

y_1 = kron(states(1,1:10:end) + vhite(1,1:10:end), ones(1,10));
y_1 = y_1(1:end-9);
y_2 = kron(states(3,1:10:end) + vhite(2,1:10:end), ones(1,10));
y_2 = y_2(1:end-9);

figure(4);
subplot(4,1,1);
plot(t, states(1,:), t, y_1);
xlabel('Time [s]');
ylabel('Position x [m]');
legend('True', 'Measured');
grid on;
%saveas(gcf, 'position_task3', 'epsc');

subplot(4,1,2);
plot(t, a, t, a_sig);
xlabel('Time [s]');
ylabel('Acceleration a [m/s^2]');
grid on;
%saveas(gcf, 'accel_task3', 'epsc');

subplot(4,1,3);
plot(t, states(3,:), t, y_2);
xlabel('Time [s]');
ylabel('Angle \theta [rad]');
grid on;
%saveas(gcf, 'orientation_task3', 'epsc');

subplot(4,1,4);
plot(t, u(2,:), t, omega_sig);
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');
grid on;
saveas(gcf, 'angular_vel_task3', 'epsc');

%% Kalman filter
x_hat = disc_dir_kalman(u, t, white, vhite, [y_1 ; y_2]);

figure(5);
subplot(5,1,1);
plot(t, states(1,:), t, x_hat(1,:));
legend('True', 'Estimated');
title('Position');
grid on;

subplot(5,1,2);
plot(t, states(2,:), t, x_hat(2,:));
title('Velocity');
grid on;

subplot(5,1,3);
plot(t, bias(1,:), t, x_hat(3,:));
title('b_1');
grid on;

subplot(5,1,4);
plot(t, states(3,:), t, x_hat(4,:));
title('Orientation \theta');
grid on;

subplot(5,1,5);
plot(t, bias(2,:), t, x_hat(5,:));
title('b_2');
grid on;