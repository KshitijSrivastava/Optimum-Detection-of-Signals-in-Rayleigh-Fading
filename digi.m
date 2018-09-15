clear all;
clc;
sig=0.01; %Sigma square
count=1;
for i=0:0.0035:35
    e_b=i/1;
    %p_b_coh(count)=1/2*(1-sqrt(e_b/(e_b+2)));
    
    pb_e1(count)=(1/2)*(1-(sqrt((e_b*(e_b*sig+2-2*sig))/((e_b+2)*(e_b*sig+2)))));
    pb_e2(count)=(1/2)*(1-sqrt((e_b)/(2+e_b)));
    pb_e3(count)=1/(2+e_b);
    pb_e4(count)=(1/2)*(1-sqrt(((1-sig)*e_b)/(1+e_b)));
    
    
    x(count)=e_b;
    count=count+1;
end

semilogy(x,pb_e1);
hold on;
grid on;
semilogy(x,pb_e2);
semilogy(x,pb_e3);
semilogy(x,pb_e4);
hold off;



% semilogy(db(x),pb_e1);
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')
% hold on;
% grid on;
% semilogy(db(x),pb_e2);
% semilogy(db(x),pb_e3);
% semilogy(db(x),pb_e4);
% hold off;
