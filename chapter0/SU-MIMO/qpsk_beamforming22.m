% qpsk BER�������
% AWGN��Rayleigh˥���ŵ�MIMOϵͳ���棨�������ζԱȰ汾��

%% ����׼��
clear all; close all; clc;
N = 2^21;
m = N/2;
maxlen = sqrt(2);
nTx = 2;
nRx = 2;
Eb_N0_dB = -2:1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'��ȴ������ڼ���...');

%% ��ʼ����
for ii = 1:length(Eb_N0_dB)
   % ����qpsk�ź�
   ip = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % ʹ�ñ���������qpsk����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(2,m) + 1j*randn(2,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nTx*nRx,m) + 1j*randn(nTx*nRx,m)); % �����ŵ�
   hEff = h.*exp(-1j*angle(h)); % ���Ͷ˲�������
   sigma=calnoisesigma(4,s,Eb_N0_dB(ii),1); % ��������
   sr = (1/sqrt(nTx))*kron(ones(nTx*nRx,1),s); 
   %�����������ߵķ������� ��Ϊһ��4*N�ľ���
   
   % ģ�⾭��AWGN�����ŵ�
   sray=hEff.*sr;
   y2 = [sray(1,:)+sray(2,:);sray(3,:)+sray(4,:)] + sigma*n; 
   y2 = sum(y2);
   
   % ���ն˾�����
   y2Hat = (maxlen).*y2./sum(hEff,1); 
   %figure; scatter(real(y1Hat),imag(y1Hat));
   
   % ���ն�Ӳ�о�
   [ip2Hat,uu2] = qpsk_demodulation(y2Hat);

   % ���������ߺͷ�����������
   nErr2(ii) = size(find(ip- ip2Hat),2);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',4,4);
end
simBer2 = nErr2/N; % ����BER

%% Alamouti
precision=2^21;
sim=zeros(1,length(Eb_N0_dB));
for ii=1:length(Eb_N0_dB)
    sim(ii)=alamouti_sim_mimo(precision,4,Eb_N0_dB(ii));
end

%% ��ͼ
close(hh);
figure;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2); hold on;
semilogy(Eb_N0_dB,sim,'r*','LineWidth',2); hold on;
axis normal tight; grid on;
legend('L=4 ������','2��2�� �������η���','2��2�� Alamouti����');
xlabel('EbN0 (dB)'); ylabel('BER');
title('QPSK��Beamforming��Alamouti���� (˫��˫��MIMO)');