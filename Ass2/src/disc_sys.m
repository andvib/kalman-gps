function x = disc_sys(x, u, t)
    h = 0.01; % [s]

    A = [[1 h 0];
         [0 1 0];
         [0 0 1]];
     
    B = [[0 0];
         [h 0];
         [0 h]];
    
    for i = 1:length(t)-1
        x(:,i+1) = A*x(:,i) + B*u(:,i);
    end
end