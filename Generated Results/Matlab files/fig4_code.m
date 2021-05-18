
beta = 1/3;
m = 1;
n = 5;
M = 5;
rho_m = [10 10 10 10 10 10 10];
rho_n = [0 5 10 15 20 25 30];
for itr = 1:length(rho_m)
    rho_m_lin = db2pow(rho_m(itr));
    rho_n_lin = db2pow(rho_n(itr));
    
    c_mn = factorial(M) / (factorial(m-1)* factorial(n-1-m) * factorial(M-n));
    
    
    out_added = 0;
    t1_sum = 0;
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
            out_added = c_p * added + out_added;
        end
        Pn_hat(itr) = 1 - (c_mn * out_added);
        disp(Pn_hat);
    else
        k1 = (1 - 2*beta) /((rho_n_lin*power(beta, 2)) - ((1 - beta)*rho_m_lin));
        for i =1: m
            l = i-1;
            c_l = nchoosek(m-1, l) * power(-1, l);
            t1_ad = exp((-(M - m  + l + 1)*k1)/ (M - m  + l + 1));
            t1_sum = t1_ad*c_l + t1_sum;
        end
        T1 = t1_sum * factorial(M) / (factorial(M -m) * factorial(m -1));
        
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
            T2_out_added = c_p * added + out_added;
        end
        T2 = T2_out_added * c_mn;
        
        Pn_hat(itr) = 1 - T1 -T2
    end
end

close all
semilogy([0 5 10 15 20 25 30], Pn_hat, 'go-', 'LineWidth', 1);
hold on;
axis([0 20 10^(-4) 1])
xlabel('User n''s transmit power in dB (rho_n)');
ylabel('Probability for OMA-MEC to outperform NOMA-MEC');
grid on
legend('m=3, n=5, beta=1/4, ana.');




