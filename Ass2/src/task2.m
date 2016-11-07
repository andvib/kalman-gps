close all;
clear all;

t = [0:0.01:20];

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
title('Acceleration a')

subplot(3,1,3);
plot(t, states(3,:));
grid on;
ylabel('[rad/s]'),
xlabel('[s]')
title('Angular velocity \theta')

saveas(gcf, 'disc_states', 'epsc');

%% Simulate discrete bias and white noise
[bias, white] = disc_bias(t);

figure(3);
plot(t, bias);
legend('b_1', 'b_2');
xlabel('[t]');
grid on;
saveas(gcf, 'disc_bias', 'epsc');

%% Simulate system with noise and bias
a_sig = states(2,:) + bias(1,:) + white(1,:);
omega_sig = states(3,:) + bias(2,:) + white(2,:);

states_sig = zeros(3,length(t));
states_sig = disc_sys(states_sig, [a_sig ; omega_sig], t);

figure(4);
plot(t, states(1,:), t, states_sig(1,:));
xlabel('Time [s]');
ylabel('Position x [m]');
legend('True', 'Measured');
grid on;
saveas(gcf, 'position_task3', 'epsc');

figure(5);
plot(t, u(1, :), t, a_sig);
xlabel('Time [s]');
ylabel('Acceleration a [m/s^2]');
legend('True', 'Measured');
grid on;
saveas(gcf, 'accel_task3', 'epsc');

figure(6);
plot(t, states(3,:), t, states_sig(3,:));
xlabel('Time [s]');
ylabel('Angle \theta [rad]');
legend('True', 'Measured');
grid on;
saveas(gcf, 'orientation_task3', 'epsc');

figure(7);
plot(t, u(2,:), t, omega_sig);
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');
legend('True', 'Measured');
grid on;
saveas(gcf, 'angular_vel_task3', 'epsc');