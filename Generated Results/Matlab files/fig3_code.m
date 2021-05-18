clc;
close all;
m = 1;
n = 5;
M = 5;

rho_n = 10:5:70;

rho_n_lin = db2pow(rho_n);
eta = 1;
rho_m_lin = rho_n_lin / eta;
c_mn = factorial(M) / (factorial(m-1)* factorial(n-1-m) * factorial(M-n));

P_n = zeros([1 length(rho_n_lin)]);
for itr = 1:length(rho_n)
    added = 0;
    for i = 1: n - m
        p = i -1;
        c_p = nchoosek(n-1-m, p) * power(-1, n - 1 - m -p);
        t1 = power(eta, m/2) * c_mn * c_p;
        t2 = power(M - m -p, (m/2) + 1);
        
        lambda = p + 1 + (M - m -p)/(eta);
        mu_m = 0;
        for j = 1:m
            l = j-1;
            c_l = nchoosek(m-1, l) * power(-1, l);
            ut1 = c_l * power(l, m);
            ut2 = factorial(m) * lambda * power(-1, m-1);
            mu_m = mu_m + ut1 + ut2;
        end
        
        
        a = (power(rho_m_lin(itr), 2) / rho_n_lin(itr)) * (M - m - p);   
        Q1 = 0; Q2 = 0;
        Q1t1 = 0;
        Q1t2 = 0;
        Q2t1 = 0;
        Q2t2 = 0;
        if mod(m,2)==1
            Q1t1 = power(pi, 0.5) * power(-1, m - 1) * factorial(m-1);
            Q1t2 = factorial((m-1)/2)*power(2, m);
            Q1 = Q1t1 / Q1t2;
            
            Q2t1 = mu_m;
            Q2t2 = double_factorial(m) * power(2, (m+1)/2) * power(a, 0.5);
            
            Q2 = Q2t1 / Q2t2;
        else
            Q1t1 = mu_m * power(pi, 0.5);
            Q1t2 = factorial(m/2) * power(2, m+1) * power(a , 0.5);
            Q1 = Q1t1  / Q1t2;
            
            Q2t1 = power(-1, m-1) * factorial(m-1);
            Q2t2 = double_factorial(m-1) * power(2, m/2);
            Q2 = Q2t1 / Q2t2;
        end
        
        added = added + ((t1/t2) * (Q1 - Q2));
        
    end
    
    P_n(itr) = added / power(rho_m_lin(itr), m/2);
    disp(P_n(itr));
end

semilogy(rho_n, P_n,'o-');
legend('m=1 n=5 eta=1');
xlabel('user n transmit power, dB');
ylabel('offloading probability');
axis([ 30 70 10^(-5) 1]);
