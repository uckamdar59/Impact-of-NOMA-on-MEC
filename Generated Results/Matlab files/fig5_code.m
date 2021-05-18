beta = 1/3;
m = 1;
n = 5;
M = 5;
eta = 5;
rho_n = [0 5 10 15 20 25];
for itr = 1:length(rho_m)
    rho_n_lin = db2pow(rho_n(itr));
    rho_m_lin = rho_n_lin / eta;
    c_mn = factorial(M) / (factorial(m-1)* factorial(n-1-m) * factorial(M-n));
    
    
    output = 0;
    temp_sum = 0;
    if (1 - beta)*rho_m_lin >= rho_n_lin * (beta.^2)
        
        for i = 1:n - m
            p = i -1;
            c_p = nchoosek(n-1-m, p) *power(-1, n - 1 - m -p);
            added = 0;
            for j = 1: m
                l = j -1;
                c_l = nchoosek(m-1, l) * power(-1, l);
                ahat = (rho_m_lin * (1 - beta) * (M  -m -p))/ (rho_n_lin *power(beta, 2)) + p + l +1;
                t1 = exp(-(M - m -p) * ((1 - 2*beta) /(rho_n_lin * power(beta, 2))));
                t2 = (M -m -p) * ahat;
                added = (c_l* t1/t2) + added;
            end
            output = c_p * added + output;
        end
        Pn_hat(itr) = 1 - (c_mn * output);
        disp(Pn_hat);
    else
        k1 = (1 - 2*beta) /((rho_n_lin*power(beta, 2)) - ((1 - beta)*rho_m_lin));
        for i =1: m
            l = i-1;
            c_l = nchoosek(m-1, l) * power(-1, l);
            t1_ad = exp((-(M - m  + l + 1)*k1)/ (M - m  + l + 1));
            temp_sum = t1_ad*c_l + temp_sum;
        end
        T1 = temp_sum * factorial(M) / (factorial(M -m) * factorial(m -1));
        
        for i = 1:n - m
            p = i -1;
            c_p = nchoosek(n-1-m, p) *power(-1, n - 1 - m -p);
            added = 0;
            for j = 1: m
                l = j -1;
                c_l = nchoosek(m-1, l) * power(-1, l);
                ahat = (rho_m_lin * (1 - beta) * (M  -m -p))/ (rho_n_lin *power(beta, 2)) + p + l +1;
                t1 = exp(-(M - m -p) * ((1 - 2*beta) /(rho_n_lin * power(beta, 2))));
                t2 = (M -m -p) * ahat;
                t3 = 1 - exp(-ahat*k1);
                T2_added = ((c_l*t1*t3)/t2) + added;
            end
            T2_out_added = c_p * added + output;
        end
        T2 = T2_out_added * c_mn;
        
        Pn_hat(itr) = 1 - T1 -T2;
    end
end

semilogy([0 5 10 15 20 25], Pn_hat, 'rx-');
axis([0 25 10^(-5) 1])
xlabel('User n''s transmit power in dB (rho_n)');
ylabel('Probability for OMA-MEC to outperform NOMA-MEC');
grid on
legend('m=1, n=5, beta=1/3');