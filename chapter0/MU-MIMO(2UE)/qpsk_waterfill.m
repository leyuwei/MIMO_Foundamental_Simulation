% qpsk BER�������
% AWGN��Rayleigh˥���ŵ�MIMOϵͳ����

%% ����׼��
clear all; close all; clc;
format long;
N = 2^18;
m = N/2;
maxlen = sqrt(2);
nTx = 2;
nRx = 2;
Eb_N0_dB = -4:1:13;
nErr2 = zeros(1,length(Eb_N0_dB));
nErr2w = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));

%% �ȹ��ʷ���
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % ʹ�ñ�����1����qam����ͼ
   s2 = qpsk_modulation(ip2); % ʹ�ñ�����2����qam����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   s2 = s2./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   snr=Eb_N0_dB(ii)+10*log10(2);
   sigma=1/(10^(snr/10));
   sigma=sqrt(sigma);
   sr = zeros(nTx,nRx,m);
   for t=1:nTx
        sr(t,1,:)=(1/sqrt(nTx))*s; 
        sr(t,2,:)=(1/sqrt(nTx))*s2; 
   end
   %�����������ߵķ������� ��Ϊһ��2*2*m�ľ���
   
   % ģ�⾭��AWGN�����ŵ�
   % CSI�ֽ�
   [precode, equalizer, diagonal, sray]=svdChannel(h,sr);
   y = sray + sigma.*n; 
   
   % ���ն˾�����
   y2Hat=zeros(nRx,nRx,m);
   for jj=1:m
       y2Hat(:,:,jj) = equalizer(:,:,jj)*y(:,:,jj);
   end
   
   % ���ն�MIMO�źż��
   tz_r_sum=zeros(1,nRx,m);
   rr=zeros(1,1,m);
   for t=1:nTx
       for r=1:nRx
           tz_r_sum(1,r,:)=tz_r_sum(1,r,:)+y2Hat(t,r,:).*conj(diagonal(t,t,:));
       end
       rr(1,1,:)=rr(1,1,:)+abs(diagonal(t,t,:)).^2;
   end
   for r=1:nRx
       tz_rr(r,:)=tz_r_sum(1,r,:)./rr(1,1,:);
   end
   recvsig1 = tz_rr(1,:);
   recvsig2 = tz_rr(2,:);
   y2Hat = recvsig1.*(sqrt(nTx));
   y2Hat = reshape(y2Hat,1,m);
   y2Hat2 = recvsig2.*(sqrt(nTx));
   y2Hat2 = reshape(y2Hat2,1,m);
   
   % ���ն�Ӳ�о�
   ip2Hat = qpsk_demodulation(y2Hat);
   ip2Hat2 = qpsk_demodulation(y2Hat2);

   % ���������ߺͷ�����������
   nErr2(ii) = size(find(ip2 - ip2Hat2),2) + size(find(ip - ip2Hat),2);
   disp(['error rate: ' num2str(nErr2(ii)/(N*2))]);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',4,4);
end
simBer2 = nErr2/(N*2); % ����BER

%% עˮ�㷨
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % ʹ�ñ�����1����qam����ͼ
   s2 = qpsk_modulation(ip2); % ʹ�ñ�����2����qam����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   s2 = s2./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   snr=Eb_N0_dB(ii)+10*log10(2);
   sigma=1/(10^(snr/10));
   sigma=sqrt(sigma);
   
   % עˮ����
   allocatedpwr=zeros(nTx,m);
   waterlv=zeros(1,m);
   ebn0=10^(Eb_N0_dB(ii)/10);
   for k=1:m
       [~,~,diagonal]=svdChannel(h(:,:,k)'*h(:,:,k));
       channelsigma=(diag(diagonal)).^2;
       bottomlevel=sigma./channelsigma;% �����������ף����ײο���
       [allocatedp,water]=waterfill2(1,bottomlevel);% ��һ��waterfill�㷨
       allocatedpwr(:,k)=allocatedp';
       waterlv(1,k)=water;
   end
   
   sr = zeros(nTx,nRx,m);
   for t=1:nTx
        sr(t,1,:)=s.*sqrt(allocatedpwr(t,:)); % ��̬���ʷ��� allocatedpwr(1,:)
        sr(t,2,:)=s2.*sqrt(allocatedpwr(t,:)); % ��̬���ʷ���
   end
   %�����������ߵķ������� ��Ϊһ��2*2*m�ľ���
   
   % ģ�⾭��AWGN�����ŵ�
   [precode, equalizer, diagonal, sray]=svdChannel(h,sr);% CSI�ֽ�
   y = sray + sigma.*n; 
   
   % ���ն˾�����
   y2Hat=zeros(nRx,nRx,m);
   for jj=1:m
       y2Hat(:,:,jj) = equalizer(:,:,jj)*y(:,:,jj);
   end
   
   % ���ն�MIMO�źż��
   tz_r_sum=zeros(1,nRx,m);
   rr=zeros(1,1,m);
   for t=1:nTx
       for r=1:nRx
           tz_r_sum(1,r,:)=tz_r_sum(1,r,:)+y2Hat(t,r,:).*conj(diagonal(t,t,:));
       end
       rr(1,1,:)=rr(1,1,:)+abs(diagonal(t,t,:)).^2;
   end
   for r=1:nRx
       tz_rr(r,:)=tz_r_sum(1,r,:)./rr(1,1,:);
   end
   recvsig1 = tz_rr(1,:);
   recvsig2 = tz_rr(2,:);
   y2Hat = recvsig1.*(sqrt(nTx));
   y2Hat = reshape(y2Hat,1,m);
   y2Hat2 = recvsig2.*(sqrt(nTx));
   y2Hat2 = reshape(y2Hat2,1,m);
   
   % ���ն�Ӳ�о�
   ip2Hat = qpsk_demodulation(y2Hat);
   ip2Hat2 = qpsk_demodulation(y2Hat2);

   % ���������ߺͷ�����������
   nErr2w(ii) = size(find(ip2 - ip2Hat2),2) + size(find(ip - ip2Hat),2);
   disp(['error rate (water-filling): ' num2str(nErr2w(ii)/(N*2))]);
end
sim = nErr2w/(N*2); % ����BER

%% ��ͼ
figure;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2); hold on;
semilogy(Eb_N0_dB,sim,'r*','LineWidth',2); hold on;
axis normal tight; grid on;
legend('4������ֵ','2��2�� �ȹ��ʷ���','2��2�� עˮ�㷨����');
xlabel('EbN0 (dB)'); ylabel('BER');
title('QPSK�ĵȹ��ʺ�עˮ������� (MIMO)');