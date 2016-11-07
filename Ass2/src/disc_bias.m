function [bias, white] = disc_bias(t)
    %% Generate white noise
    w_1 = 0.2 * wgn(1, length(t), 1);
    w_2 = 0.8 * wgn(1, length(t), 1);
    w_3 = 0.2 * wgn(1, length(t), 1);
    w_4 = 0.8 * wgn(1, length(t), 1);
    
    %% Simulate system
    x = zeros(2, length(t));
    
    T_1 = 10;
    T_2 = 10;
    h = 0.01;
    
    A = [[(1-(h/T_1)) 0];
         [0 (1-(h/T_2))]];
    E = [[h 0] ; [0 h]];
    
    for i = 1:length(t)-1
        x(:,i+1) = A*x(:,i) + E*[w_2(i) ; w_4(i)];
    end
    
    bias = x;
    white = [w_1; w_3];
end