K = 4;
Pe_PSK = zeros(1, 11);
Pe_QAM = zeros(1, 11);
Pe_PPM = zeros(1, 11);
for n_r = 1:11
    N_s = (n_r-1) * log2(K);
    % 4-PSK
    delta_PSK = sqrt(N_s);
    lambda = zeros(1, K);
    % find lambdas 
    for i = 1 : K
        for j = 1 : K
          lambda(i) = lambda(i) + gram_psk(0, (j- 1), delta_PSK, K) * exp(-1i* (i-1) * (j-1) * 2 * pi/K);
        end
    end
    
    % sum over desquared lambdas
    sum_PSK = 0;
    for i = 1:K
        sum_PSK = sum_PSK + sqrt(lambda(i));
    end

    Pe_PSK(n_r) = real(1 - ((sum_PSK)^2)/K^2);

    % 4-QAM
    delta_QAM = sqrt((3/2 * N_s)/(K - 1));
    L = sqrt(K);
    G = zeros(K, K);
    for i =  1 : L
        for j = 1 : L
            u = -(L-1)+2*(i-1);
            v = -(L-1)+2*(j-1); 
            for p = 1 : L
                for q = 1 : L
                    u_pr = -(L-1)+2*(p-1);
                    v_pr = -(L-1)+2*(q-1);
                    alpha = u + 1i * v;
                    beta = u_pr + 1i * v_pr;
                    G(2*(i-1) + j, 2 * (p-1) + q) = exp(- (delta_QAM^2)/2 *(...
                      abs(alpha)^2 + abs(beta)^2 - 2 * conj(alpha) * beta));
                end
            end
        end
    end
   
    half_G = sqrtm(G);
    sum_QAM = 0;
    for i = 1 : K 
        sum_QAM = sum_QAM + abs(half_G(i, i))^2;
    end
   
    Pe_QAM(n_r) = 1 - sum_QAM/K;


    % 4-PPM
    X_2 = exp(-N_s);
    Pe_PPM(n_r) = 1 - ((sqrt(1 + (K -1 ) * X_2) + (K - 1) * sqrt(1 - X_2) )^2)/(K^2);
    
end
plot(0:10,log10(Pe_PSK))

hold on 
plot(0:10, log10(Pe_QAM))
plot(0:10, log10(Pe_PPM))

hold off















function  gram_pq = gram_psk(p, q, del, Ka)
    gram_pq = exp(-(del^2)*(1 - exp(1i* (q-p) * 2 * pi/Ka)));
end