% QPSK BER�������
% AWGN��Rayleigh˥���ŵ�MISOϵͳ���棨�������ζԱȰ汾��

%% ����׼��
clear all; close all; clc;
N = 2^16;
m = N/2;
maxlen = sqrt(2);
nTx = 2;
Eb_N0_dB = -2:1:20;
nErr1 = zeros(1,length(Eb_N0_dB));
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'��ȴ������ڼ���...');

%% ��ʼ����
for ii = 1:length(Eb_N0_dB)
   % ����qpsk�ź�
   ip = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % ʹ�ñ���������qpsk����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % �����ŵ�
   hEff = h.*exp(-1j*angle(h)); % ���Ͷ˲�������
   sigma=calNoisePower(4,s,Eb_N0_dB(ii)); % ��������
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %�����������ߵķ������� ��Ϊһ��2*N�ľ��� Ҳ����˵��һ����ά����
   
   % ģ�⾭��AWGN�����ŵ�
   y1 = sum(h.*sr,1) + sigma*n; 
   y2 = sum(hEff.*sr,1) + sigma*n; 

   % ���ն˾�����
   y1Hat = maxlen.*y1./sum(h,1); 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   % figure; scatter(real(y1Hat),imag(y1Hat));
   
   % ���ն�Ӳ�о�
   [ip1Hat,uu] = qpsk_demodulation(y1Hat);
   [ip2Hat,uu2] = qpsk_demodulation(y2Hat);

   % ���������ߺͷ�����������
   nErr1(ii) = size(find(ip- ip1Hat),2);
   nErr2(ii) = size(find(ip- ip2Hat),2);
   theoryBer_nRx2(ii) = berfading(Eb_N0_dB(ii),'psk',4,1);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',4,2);
end
simBer1 = nErr1/N; % ����BER ������Beamforming��
simBer2 = nErr2/N; % ����BER ����Beamforming��

%% ��ͼ
close(hh);
figure;
semilogy(Eb_N0_dB,theoryBer_nRx2,'b*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer1,'r*','LineWidth',2); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2);
axis normal tight; grid on;
legend('L=1 ������','L=2 ������','2��1�� �޲������η���','2��1�� �������η���');
xlabel('EbN0 (dB)'); ylabel('BER');
title('QPSK��AWGN����˥���ŵ��µ�MISOϵͳ���� (��������)');