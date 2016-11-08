function x_hat = disc_dir_kalman(u, t, w, v, y)
    %T_1 = 0.01;
    %T_2 = 0.01;
    global T_1 T_2;
    
    %% Create matrices
    h = 0.01;
    Phi = [[1 h 0 0 0];
           [0 1 h 0 0];
           [0 0 (1-(h/T_1)) 0 0];
           [0 0 0 1 h];
           [0 0 0 0 (1-(h/T_2))]];
       
    delta = [[0 0];
             [h 0];
             [0 0];
             [0 h];
             [0 0]];
         
    H = [[1 0 0 0 0];
         [0 0 0 1 0]];
     
    Gamma = [[0 0 0 0];
             [1 0 0 0];
             [0 1 0 0];
             [0 0 1 0];
             [0 0 0 1]];
    
    %% Create design matrices
    Q = h * diag([var(w(1,:)), var(w(2,:)), var(w(3,:)), var(w(4,:))]);
    R = (1/h) * diag([var(v(1,:)), var(v(2,:))]);

    %% Initialize vectors
    x_bar = zeros(5, length(t));
    x_bar(:,1) = [0; 0; 0; 0; 0];
    
    x_hat = zeros(5, length(t));
    
    P_bar = eye(5);
    P_bar(1,1) = 0;
    %% Main loop
    for i = 1:length(t)-1
        % Calculate Kalman gain matrix
        K = P_bar*H'*inv(H*P_bar*H' + R);
        
        % State estimation update
        x_hat(:,i) = x_bar(:,i) + K*(y(:,i) - H*x_bar(:,i));
        
        % Error covariance update
        P_hat = (eye(5)-K*H) * P_bar * (eye(5)-K*H)' + K*R*K';
        
        % State estimation propagation
        x_bar(:,i+1) = Phi*x_hat(:,i) + delta*u(:,i);
        
        % Error covariance propagation
        P_bar = Phi*P_hat*Phi' + Gamma*Q*Gamma';
    end
end