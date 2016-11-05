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