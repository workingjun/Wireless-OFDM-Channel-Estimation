function [out] = Compare_fftshift(Before_fftshift, After_fftshift)
%UNTITLED 이 함수의 요약 설명 위치
%   자세한 설명 위치
figure(1004);

plot(abs(Before_fftshift).^2, 'bo-');
hold on;
plot(abs(After_fftshift).^2, 'rx-');
grid on;
legend('Before_fftshift', 'After_fftshift', 'Location','southwest');
% xlabel('Real of Rx Signal', 'Fontsize', 12);
% ylabel('Imag of Rx Signal', 'Fontsize', 12);
% title('Constellation for BPSK', 'Fontsize', 14);
%axis([-3 3 -3 3]);
pause();

out = 1;

