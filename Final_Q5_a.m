% a 
N_set = [0, 0.05, 0.1, 0.2];
Pe = zeros(4, 17);
num = 40;
counter = 1;
for N = N_set
    for Ns = 0.5:0.5:8
        if N == 0 
            Pe(counter, 2*Ns + 1) = (1 - sqrt(1 - exp(-2*Ns)))/2;
        else
            alpha = sqrt(2*Ns);
            R_0 = zeros(num, num);
            for i = 1: num
                R_0(i, i) = ((N/(N+1))^(i-1))/(N + 1);
            end
            
            R_alpha = zeros(num, num);
            for i = 1:num
                for j = 1:num 
                    n = i - 1;
                    m = j - 1;
                    R_alpha(j,i) = (N^n/(N+1)^(n+1)) * sqrt(factorial(m)/factorial(n)) * (conj(alpha)/N)^(n-m) * exp(-abs(alpha)^2/(N+1)) * laguerreL(m,n-m,- abs(alpha)^2/(N*(N+1)));
                end
            end

            Decision = (R_alpha - R_0)/2;
            D = eig(Decision);
            posSum = 0;
            for i = 1:numel(D)
                if D(i) > 0 
                    posSum = posSum + D(i);
                end
            end
            
            Pe(counter, 2 * Ns + 1) = 1/2 - posSum;
        end
    end
    counter = counter + 1;
end
Ns_set = 0:0.5:8;
plot( Ns_set, log10(Pe(1, :)))
hold on 
plot( Ns_set, log10(Pe(2, :)))
plot( Ns_set, log10(Pe(3, :)))
plot( Ns_set, log10(Pe(4, :)))
hold off