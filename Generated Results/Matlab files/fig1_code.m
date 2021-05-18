close all;
clear all;
clc;
M=5;
m=4;
n=5;
R_n1=[0:1:30];
R_n=10.^(R_n1/10);
R_m1=5;
R_m=10.^(R_m1/10);


C_mn=(factorial(M))./(factorial(m-1).*factorial(n-1-m).*factorial(M-n))

%------

% for ii=1:length(R_n1)
% for p=0:(n-1-m)
%     cpp= (nchoosek(n-1-m,p).*((-1)^(n-1-m-p)));
%     cp=cpp./(M-n-p);
%     
%     for l=0:(m-1)
%         a=((R_m.^2)./ii).*(M-m-p)
%         b=p+l+1+(M-m-p).*(R_m./ii)
%         cll=nchoosek(m-1,l).*((-1).^l);
%         cl=cll.*(exp(b.^2/(4.*a)));
%         
%         Q_function=qfunc((max(0,ii-R_m)./R_m.^2)+b./(2.*sqrt(a)));
%         Q_function_1=cl*Q_function*(sqrt(pi)./(2.*sqrt(a)))
% %        temp1=(factorial(M))./(factorial(m-1).*factorial(M-m));
% %       temp2=Q_function_1-temp1;
%         
%     end
%     temp3=cp*Q_function_1+1;
%     temp1=(factorial(M))./(factorial(m-1).*factorial(M-m));
%     
%     for l=0:m-1
%         temp4=exp(-(M-m+l+1).*(max(0,ii-R_m)./R_m.^2));
%         temp5=temp4./(M-m+l+1);
%         temp6=temp5.*cll; 
%     end
%     
%     temp7=temp6.*temp1;
%     temp=temp3-temp7;
%     tem(ii)=temp*C_mn
% end
% end
%---------------




for ii=1:length(R_n1)
    syms v
        f2=symsum(((nchoosek(m-1,v).*((-1).^v)).*(exp((-(M-m+v+1)).*(max(0,R_n-R_m)./R_m.^2))./((M-m+v+1)))),v,0,m-1);
    F2=(factorial(M))./(factorial(m-1).*factorial(M-m)).*f2;
    
    F3=1-F2
    
 
end


close all
figure
semilogy(R_n1,F3,'rd-','LineWidth',1);
hold on
xlabel('user n transmit power dB');
ylabel('offloading probability');
legend('m=4 n=5');

grid on