function [out] = H_Frequency_domain(h, H)
%UNTITLED �� �Լ��� ��� ���� ��ġ
%   �ڼ��� ���� ��ġ
%
% h :: �ð��� ä�� ��
% H :: ���ļ� �� ä�� ��

N = length(H);
t_index = 1:N;
f_index = -N/2:N/2-1;
h_t = zeros(1,N);
h_t(1:length(h)) = h; 


figure(1000);

subplot(2,1,1);
bar(t_index, abs(h_t).^2);%% Simulation !!
hold on;
% semilogy(SNR_dB, Sim_SER, 'ro-');%% Simulation !!
% semilogy(SNR_dB, Analysis_BER,'k*-');%%  ���ɺм� ���
% semilogy(SNR_dB, Analysis_SER,'kh-');%%  ���ɺм� ���
% legend('BER, Simulation', 'SER, Simulation' , 'BER, Analysis','SER, Analysis', 'Location','southwest');
grid on;
xlabel('Time', 'Fontsize', 12);
ylabel('|h|^2', 'Fontsize', 12);
%title('MPSK (AWGN Channel) [by KKB]', 'Fontsize', 14);
xlim([0 64]);
%ylim([0 1]);

subplot(2,1,2); 
semilogy(f_index, abs(H).^2);%% Simulation !!
hold on;
% semilogy(SNR_dB, Sim_SER, 'ro-');%% Simulation !!
% semilogy(SNR_dB, Analysis_BER,'k*-');%%  ���ɺм� ���
% semilogy(SNR_dB, Analysis_SER,'kh-');%%  ���ɺм� ���
% legend('BER, Simulation', 'SER, Simulation' , 'BER, Analysis','SER, Analysis', 'Location','southwest');
grid on;
xlabel('Frequency', 'Fontsize', 12);
ylabel('|H|^2 (dB)', 'Fontsize', 12);
%title('MPSK (AWGN Channel) [by KKB]', 'Fontsize', 14);
xlim([min(f_index) max(f_index)]);
ylim([10^-2 10]);
pause();

out = 1;
end

