N_gamma = 5;
N = 0.1;
epsilon = 10^-5;
nu = 10^-5;
gamma = sqrt(N_gamma);


for n_e = 25:-1:1
    R = zeros(n_e, n_e);
 
    for i = 1:n_e 
        for j = 1:n_e
            m = j - 1;
            n = i - 1;
            R(j,i) = (N^n)/((N+1)^(n+1)) * sqrt(factorial(m)/factorial(n)) * (conj(gamma)/N)^(n-m) * exp(-(abs(gamma)^2)/(N+1))* laguerreL(m, n-m, - abs(gamma)^2/(N*(N+1)));
        end
    end

    if trace(R) < (1 - epsilon)
        n_epsilon = n_e + 1 ;
        break
    end
end 


R = zeros(n_epsilon, n_epsilon);
 
for i = 1:n_epsilon 
    for j = 1:n_epsilon
        m = j - 1;
        n = i - 1;
        R(j,i) = (N^n)/((N+1)^(n+1)) * sqrt(factorial(m)/factorial(n)) * (conj(gamma)/N)^(n-m) * exp(-(abs(gamma)^2)/(N+1))* laguerreL(m, n-m, - abs(gamma)^2/(N*(N+1)));
    end
end


[Z, D] = eig(R);

for n_n = 6:-1:1
    Z_h = Z(:, 1:n_n);
    D_h = D(1:n_n, 1:n_n);
    beta = Z_h * sqrt(D_h);
    if immse(R, beta * ctranspose(beta)) > nu
        n_nu = n_n + 1;
        break;
    end
end

n_epsilon
n_nu
