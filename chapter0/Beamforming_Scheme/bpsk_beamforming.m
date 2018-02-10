% BPSK BER�������
% AWGN��Rayleigh˥���ŵ�MISOϵͳ���棨�������ζԱȰ汾��

%% ����׼��
clear all; close all; clc;
N = 2^16;
nTx = 2;
Eb_N0_dB = -3:1:20;
nErr1 = zeros(1,length(Eb_N0_dB));
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));

%% ��ʼ����
for ii = 1:length(Eb_N0_dB)
   ip = rand(1,N)>0.5;
   s = 2*ip-1;
   n = 1/sqrt(2)*(randn(1,N) + 1j*randn(1,N)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx,N) + 1j*randn(nTx,N)); % �����ŵ�
   sigma=calNoisePower(2,s,Eb_N0_dB(ii)); % ��������
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %�����������ߵķ������� ��Ϊһ��2*N�ľ��� Ҳ����˵��һ����ά����
   
   % ģ�⾭��AWGN�����ŵ�
   hEff = h.*exp(-1j*angle(h));
   y1 = sum(h.*sr,1) + sigma*n; 
   y2 = sum(hEff.*sr,1) + sigma*n; 

   % ���ն˾�����
   y1Hat = y1./sum(h,1); 
   y2Hat = y2./sum(hEff,1); 

   % ���ն�Ӳ�о�
   ip1Hat = real(y1Hat)>0;
   ip2Hat = real(y2Hat)>0;

   % ���������ߺͷ�����������
   nErr1(ii) = size(find([ip- ip1Hat]),2);
   nErr2(ii) = size(find([ip- ip2Hat]),2);
   theoryBer_nRx2(ii) = berfading(Eb_N0_dB(ii),'psk',2,1);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',2,2);
end
simBer1 = nErr1/N; % ����BER ������Beamforming��
simBer2 = nErr2/N; % ����BER ����Beamforming��

%% ��ͼ
figure;
semilogy(Eb_N0_dB,theoryBer_nRx2,'b*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer1,'r*','LineWidth',2); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2);
axis normal tight; grid on;
legend('L=1 ������','L=2 ������','2��1�� �޲������η���','2��1�� �������η���');
xlabel('EbN0 (dB)'); ylabel('BER');
title('BPSK��AWGN����˥���ŵ��µ�MISOϵͳ���� (��������)');