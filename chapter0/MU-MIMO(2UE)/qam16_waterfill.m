% 16QAM BER�������
% AWGN��Rayleigh˥���ŵ�MIMOϵͳ����

%% ����׼��
clear all; close all; clc;
format long;
N = 2^15;
m = N/4;
maxlen = sqrt(10);
nTx = 2;
nRx = 2;
Eb_N0_dB = -4:1:13;
% Eb_N0_dB = 1000:-1:999;
nErr2 = zeros(1,length(Eb_N0_dB));
nErr2w = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));

%% �ȹ��ʷ���
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ�����1����qam����ͼ
   s2 = qam16_modulation(ip2); % ʹ�ñ�����2����qam����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   s2 = s2./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   snr=Eb_N0_dB(ii)+10*log10(4);
   sigma=1/(10^(snr/10));
   sigma=sqrt(sigma);
   sr = zeros(nTx,nRx,m);
   ss = [s;s2];
   for r=1:nRx
       for t=1:nTx
            sr(t,r,:)=(1/sqrt(nTx))*ss(t,:); 
       end
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
   tz_r_sum=zeros(nTx,1,m);
   rr=zeros(nTx,1,m);
   for t=1:nTx
       for r=1:nRx
           tz_r_sum(t,1,:)=tz_r_sum(t,1,:)+y2Hat(t,r,:).*conj(diagonal(t,t,:));
       end
       rr(t,1,:)=rr(t,1,:)+2.*abs(diagonal(t,t,:)).^2;
   end
   for t=1:nTx
       for r=1:nRx
           tz_rr(t,:)=tz_r_sum(t,1,:)./rr(t,1,:);
       end
   end
   recvsig1 = tz_rr(1,:);
   recvsig2 = tz_rr(2,:);
   y2Hat = recvsig1.*(sqrt(nTx));
   y2Hat = reshape(y2Hat,1,m);
   y2Hat2 = recvsig2.*(sqrt(nTx));
   y2Hat2 = reshape(y2Hat2,1,m);
   
   % ���ն�Ӳ�о�
   ip2Hat = qam16_demodulation(y2Hat,maxlen);
   ip2Hat2 = qam16_demodulation(y2Hat2,maxlen);

   % ���������ߺͷ�����������
   nErr2(ii) = size(find(ip2 - ip2Hat2),2) + size(find(ip - ip2Hat),2);
   disp(['error rate: ' num2str(nErr2(ii)/(N*2))]);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'qam',16,4);
end
simBer2 = nErr2/(N*2); % ����BER

%% עˮ�㷨
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ�����1����qam����ͼ
   s2 = qam16_modulation(ip2); % ʹ�ñ�����2����qam����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   s2 = s2./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   snr=Eb_N0_dB(ii)+10*log10(4);
   sigma=1/(10^(snr/10));
   sigma=sqrt(sigma);
   
   % עˮ����
   allocatedpwr=zeros(nTx,m);
   waterlv=zeros(1,m);
   for k=1:m
       [~,~,diagonal]=svdChannel(h(:,:,k)'*h(:,:,k));
       channelsigma=(diag(diagonal)).^2;
       bottomlevel=sigma./channelsigma;% �����������ף����ײο���
       [allocatedp,water]=waterfill2(1,bottomlevel);% һ��waterfill�㷨
       allocatedpwr(:,k)=allocatedp';
       waterlv(1,k)=water;
   end
   
   sr = zeros(nTx,nRx,m);
   ss = [s;s2];
   for r=1:nRx
       for t=1:nTx
            sr(t,r,:)=ss(r,:).*sqrt(allocatedpwr(t,:)); % ��̬���ʷ���
       end
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
   ip2Hat = qam16_demodulation(y2Hat,maxlen);
   ip2Hat2 = qam16_demodulation(y2Hat2,maxlen);

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
title('16QAM�ĵȹ��ʺ�עˮ������� (MIMO)');