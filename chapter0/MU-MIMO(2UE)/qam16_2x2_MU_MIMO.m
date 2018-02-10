% 16QAM
% 2x�������� 2x�������û� MU-MIMO DL BER�������
% �����Beamforming��ʹ��ZF-BD�㷨

%% ����׼��
clear all; close all; clc;
N = 2^17;
m = N / 4;
maxlen = sqrt(10);
nTx = 2;
nRx = 2;
Eb_N0_dB = -2:1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'��ȴ������ڷ�����...');

%% ��ʼ����
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ���������QAM����ͼ
   s2 = qam16_modulation(ip2);
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   s2 = s2./maxlen;
   n = 1/sqrt(2)*(randn(nRx,1,m) + 1j*randn(nRx,1,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   sigma = calnoisesigma(16,s,Eb_N0_dB(ii),1); % ��������
   sr = zeros(nTx,1,m);
   transig = [s;s2];
   for t = 1:nTx
       sr(t,1,:) = (1/sqrt(nTx)) * transig(t,:); 
   end
   %�����������ߵķ������� ��Ϊһ��2*1*m�ľ���
   
   % ģ�⾭��AWGN�����ŵ�
   % CSI�ֽ�
   [precode, equalizer, diagonal, sray]=svdChannel(h,sr);
   y = sray + sigma.*n; % sray = H * precode_mat * signal
   
   % ���ն˾���
   y2Hat = zeros(nRx,1,m);
   for jj = 1:m
       y2Hat(:,:,jj) = equalizer(:,:,jj) * y(:,:,jj);
   end
   
   % ���ն��źż��
   outputsig = zeros(nRx,1,m);
   for r = 1:nRx
       outputsig(r,1,:) = y2Hat(r,1,:).*conj(diagonal(r,r,:))./(abs(diagonal(r,r,:))).^2;
   end
   y2Hat = outputsig * sqrt(nTx) * maxlen;
   
   % ���ն�Ӳ�о�
   ip2Hat1 = qam16_demodulation(y2Hat(1,:),maxlen);
   ip2Hat2 = qam16_demodulation(y2Hat(2,:),maxlen);

   % ���������ߺͷ�����������
   nErr2(ii) = size(find(ip - ip2Hat1),2) + size(find(ip2 - ip2Hat2),2);
   disp(['bit error rate: ' num2str(nErr2(ii)/(nTx*N))]);
end
simBer2 = nErr2 / (2*N); % ����BER

%% ��ͼ
close(hh);
figure;
semilogy(Eb_N0_dB,simBer2,'k*-','LineWidth',1); hold on;
axis normal; grid on;
legend('��������');
xlabel('EbN0 (dB)'); ylabel('BER');
title('16QAM 2X2 MU-MIMO������');