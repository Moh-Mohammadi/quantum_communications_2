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

cvx_begin SDP 
variable X(dim, dim) hermitian
minimize(trace(X))
subject to 
X > rho0; X > rho1; X > rho2; X > rho3;...
X > rho4; X > rho5; X > rho6; X > rho7;...
X > rho8; X > rho9; X > rho10; X > rho11;...
X > rho12; X > rho13; X > rho14; X>rho15;
cvx_end
copt = cvx_optval
t = (copt); t = trace(X)
Pe = 1.0 - t





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