dim = 30;
N = 0.2;
K = 16;
Ns = 0.5;

delta = sqrt((3/2)*(Ns/(K - 1)));

rho0 = thermal_rho(dim, delta*(-3 + 1i * -3 ),N);
rho1 = thermal_rho(dim, delta*(-3 + 1i * -1 ),N);
rho2 = thermal_rho(dim, delta*(-3 + 1i * 1 ),N);
rho3 = thermal_rho(dim, delta*(-3 + 1i * 3 ),N);
rho4 = thermal_rho(dim, delta*(-1 + 1i * -3 ),N);
rho5 = thermal_rho(dim, delta*(-1 + 1i * -1 ),N);
rho6 = thermal_rho(dim, delta*(-1 + 1i * 1 ),N);
rho7 = thermal_rho(dim, delta*(-1 + 1i * 3 ),N);
rho8 = thermal_rho(dim, delta*(1 + 1i * -3 ),N);
rho9 = thermal_rho(dim, delta*(1 + 1i * -1 ),N);
rho10 = thermal_rho(dim, delta*(1 + 1i * 1 ),N);
rho11 = thermal_rho(dim, delta*(1 + 1i * 3 ),N);
rho12 = thermal_rho(dim, delta*(3 + 1i * -3 ),N);
rho13 = thermal_rho(dim, delta*(3 + 1i * -1 ),N);
rho14 = thermal_rho(dim, delta*(3 + 1i * 1 ),N);
rho15 = thermal_rho(dim, delta*(3 + 1i * 3 ),N);

[Z0, D0] = eig(rho0);
beta0 = Z0 * sqrt(D0);

[Z1, D1] = eig(rho1);
beta1 = Z1 * sqrt(D1);

[Z2, D2] = eig(rho2);
beta2 = Z2 * sqrt(D2);

[Z3, D3] = eig(rho3);
beta3 = Z3 * sqrt(D3);

[Z4, D4] = eig(rho4);
beta4 = Z4 * sqrt(D4);

[Z5, D5] = eig(rho5);
beta5 = Z5 * sqrt(D5);

[Z6, D6] = eig(rho6);
beta6 = Z6 * sqrt(D6);

[Z7, D7] = eig(rho7);
beta7 = Z7 * sqrt(D7);

[Z8, D8] = eig(rho8);
beta8 = Z8 * sqrt(D8);

[Z9, D9] = eig(rho9);
beta9 = Z9 * sqrt(D9);

[Z10, D10] = eig(rho10);
beta10 = Z10 * sqrt(D10);

[Z11, D11] = eig(rho11);
beta11 = Z11 * sqrt(D11);

[Z12, D12] = eig(rho12);
beta12 = Z12 * sqrt(D12);

[Z13, D13] = eig(rho13);
beta13 = Z13 * sqrt(D13);

[Z14, D14] = eig(rho14);
beta14 = Z14 * sqrt(D14);

[Z15, D15] = eig(rho15);
beta15 = Z15 * sqrt(D15);

Gamma = [beta0 beta1 beta2 beta3 beta4 beta5 beta6 beta7 beta8 beta9 beta10 ...
    beta11 beta12 beta13 beta14 beta15];

G = Gamma' * Gamma;
[Z, D] = eig(G);
sqrt_inv_D = D^(-1/2);
sqrt_inv_G = Z * sqrt_inv_D * Z';

Pc = 1/K * trace(sqrt_inv_G);
Pe = 1 - Pc

function dens_op = thermal_rho(num, gamma, N)
    dens_op = zeros(num, num);
    if N ~= 0
        for i = 1:num 
            for j = 1:num
                m = j - 1;
                n = i - 1;
                dens_op(j,i) = (N^n)/((N+1)^(n+1)) * sqrt(factorial(m)/factorial(n)) * (conj(gamma)/N)^(n-m) * exp(-(abs(gamma)^2)/(N+1))* laguerreL(m, n-m, - abs(gamma)^2/(N*(N+1)));
            end
        end
    else
        for i = 1:num
            dens_op(i,i) = N^(i-1)/N^i;    
        end

    end
end