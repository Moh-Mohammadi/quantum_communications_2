N = 0;
Pe = zeros(1, 15);
for Ns = 0.5:0.5:8    
    if N == 0 
        X_2 = exp(-Ns);
        Pe(2 * Ns) = 1 -(( sqrt(1 + 3 * X_2) + 3 * sqrt(1- X_2))^2)/16;
    else
        num = 30;
        delta = sqrt(Ns);
        
        % Define R_0 and R_1
        R_1 = zeros(num, num);
        for i = 1:num 
            for j = 1:num
                m = j - 1;
                n = i - 1;
                R_1(j,i) = (N^n)/((N+1)^(n+1)) * sqrt(factorial(m)/factorial(n)) * (conj(delta)/N)^(n-m) * exp(-(abs(delta)^2)/(N+1))* laguerreL(m, n-m, - abs(delta)^2/(N*(N+1)));
            end
        end
        R_0 = zeros(num, num);
        for i = 1:num 
                n = i - 1;
                R_0(i,i) = (N^n)/((N+1)^(n+1));
        end
        
        % Factorization of R_0 and R_1 
        h = 4;
        
        Z_0 = eye(num);
        Z_h_0 = Z_0(:, 1:h);
        D_h_0 = R_0(1:h, 1:h);
        beta_0 = Z_h_0 * sqrt(D_h_0);
        
        [Z_1, D_1] = eig(R_1);
        Z_h_1 = Z_1(:, 1:h);
        D_h_1 = D_1(1:h, 1:h);
        beta_1 = Z_h_1 * sqrt(D_h_1);
        
        beta0 = kron(beta_0,kron(beta_0,kron(beta_0, beta_1)));
        beta1 = kron(beta_0,kron(beta_0,kron(beta_1, beta_0)));
        beta2 = kron(beta_0,kron(beta_1,kron(beta_0, beta_0)));
        beta3 = kron(beta_1,kron(beta_0,kron(beta_0, beta_0)));
        
        half_d_sum = zeros(h^4, h^4);
        Wk = exp(1i * pi/2);
        for i = 1:4
            d = beta0' * beta0 + beta0.' * beta1 * Wk^(-(i-1))+ ...
            beta0.' * beta2 * Wk^(-2*(i-1)) + beta0.' * beta3 * Wk^(-3*(i-1));
            half_d_sum = half_d_sum + sqrtm(d);
        end
    
        
        Pc = trace((half_d_sum^2)/16);
        Pe(2 * Ns) = 1 - Pc;
    end

end

plot(0.5:0.5:8,Pe)

