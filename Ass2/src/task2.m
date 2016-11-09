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

y_1 = kron(states(1,1:10:end) + vhite(1,1:10:end), ones(1,10));
y_1 = y_1(1:end-9);
y_2 = kron(states(3,1:10:end) + vhite(2,1:10:end), ones(1,10));
y_2 = y_2(1:end-9);

figure(4);
subplot(2,1,1);
plot(t, states(1,:), t, y_1);
ylim([-5 50]);
xlabel('Time [s]');
ylabel('[m]');
title('Position x');
legend('True', 'Measured');
grid on;

subplot(2,1,2);
p = plot(t, a, t, a_sig);
xlabel('Time [s]');
ylabel('[m/s^2]');
title('Acceleration a');
grid on;

saveas(gcf, 'acc_pos_task3', 'epsc');

figure(5);
subplot(2,1,1);
plot(t, states(3,:), t, y_2);
xlabel('Time [s]');
ylabel('[rad]');
title('Angle \theta');
legend('True', 'Measured');
grid on;

subplot(2,1,2);
plot(t, u(2,:), t, omega_sig);
xlabel('Time [s]');
ylabel('[rad/s]');
title('Angular velocity \omega');
grid on;
saveas(gcf, 'angle_task3', 'epsc');

%% Direct Kalman filter
x_hat = disc_dir_kalman(u, t, white, vhite, [y_1 ; y_2]);

figure(6);
subplot(3,1,1);
plot(t, states(1,:), t, x_hat(1,:));
ylim([0 50]);
legend('True', 'Estimated', 'Location', 'NW');
title('Position');
grid on;

subplot(3,1,2);
plot(t, states(2,:), t, x_hat(2,:));
ylim([-1 2]);
title('Velocity');
grid on;

subplot(3,1,3);
plot(t, states(3,:), t, x_hat(4,:));
ylim([-1 2]);
title('Orientation \theta');
grid on;

saveas(gcf, 'kalman_task4', 'epsc');

figure(7);
subplot(2,1,1);
plot(t, bias(1,:), t, x_hat(3,:));
legend('True', 'Estimated');
title('b_1');
grid on;

subplot(2,1,2);
plot(t, bias(2,:), t, x_hat(5,:));
title('b_2');
grid on;

saveas(gcf, 'kalman_bias_task4', 'epsc');

%% Indirect Kalman filter
x_measured = zeros(3,length(t));
ins = disc_sys(x_measured, u, t);

delta_x_hat = indir_kalman(u, t, white, vhite, [(y_1 - ins(1,:)) ; (y_2 - ins(3,:))]);

figure(8);
subplot(3,1,1);
plot(t, states(1,:), t, y_1-delta_x_hat(1,:));
ylim([0 50]);
legend('True', 'Estimated', 'Location', 'NW');
title('Position');
grid on;

subplot(3,1,2);
plot(t, states(2,:), t, ins(2,:)-x_hat(2,:));
ylim([-1 2]);
title('Velocity');
grid on;

subplot(3,1,3);
plot(t, states(3,:), t, ins(3,:)-x_hat(4,:));
ylim([-1 2]);
title('Orientation \theta');
grid on;

saveas(gcf, 'kalman_task5', 'epsc');

figure(9);
subplot(2,1,1);
plot(t, bias(1,:), t, delta_x_hat(3,:));
legend('True', 'Estimated');
title('b_1');
grid on;

subplot(2,1,2);
plot(t, bias(2,:), t, delta_x_hat(5,:));
title('b_2');
grid on;

saveas(gcf, 'kalman_bias_task5', 'epsc');