% 16QAM BER�������
% AWGN��Rayleigh˥���ŵ�MIMOϵͳ���棨�������ζԱȰ汾��

%% ����׼��
clear all; close all; clc;
N = 2^16;
m = N/4;
maxlen = sqrt(10);
nTx = 2;
nRx = 2;
Eb_N0_dB = -2:1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));
%hh = waitbar(1,'��ȴ������ڼ���...');

%% ��ʼ����
for ii = 1:length(Eb_N0_dB)
   % ����qam�ź�
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % ʹ�ñ���������qam����ͼ
   s = s./maxlen; % ������һ����ֻ�Ķ���������
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % ��˹������
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % �����ŵ�
   sigma=calnoisesigma(16,s,Eb_N0_dB(ii),1); % ��������
   sr = zeros(nTx,nRx,m);
   for r=1:nRx
       for t=1:nTx
           sr(t,r,:)=(1/sqrt(nTx))*s; 
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
   tz_r_sum=zeros(1,nRx,m);
   rr=zeros(1,1,m);
   for t=1:nTx
       for r=1:nRx
           tz_r_sum(1,r,:)=tz_r_sum(1,r,:)+y2Hat(t,r,:).*conj(diagonal(t,t,:));
       end
       rr(1,1,:)=rr(1,1,:)+(abs(diagonal(t,t,:))).^2;
   end
   tz_rr=zeros(nRx,m);
   for r=1:nRx
       tz_rr(r,:)=tz_r_sum(1,r,:)./rr(1,1,:);
   end
%     for t=1:m
%         y2Hat(:,:,t)=y2Hat(:,:,t)./diagonal(:,:,t);
%     end
    
%    recvsig=[y2Hat(1,1,:)+y2Hat(2,1,:) y2Hat(1,2,:)+y2Hat(2,2,:)];
%    recvsig=sum(recvsig,2);
   recvsig = tz_rr(1,:)+tz_rr(2,:);
   y2Hat2 = recvsig./(sqrt(nTx));
   y2Hat2 = reshape(y2Hat2,1,m);
   
   % ���ն�Ӳ�о�
   ip2Hat = qam16_demodulation(y2Hat2,maxlen);

   % ���������ߺͷ�����������
   nErr2(ii) = size(find(ip- ip2Hat),2);
   disp(['error rate: ' num2str(nErr2(ii)/N)]);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'qam',16,4);
end
simBer2 = nErr2/N; % ����BER

%% Alamouti
precision=2^19;
sim=zeros(1,length(Eb_N0_dB));
for ii=1:length(Eb_N0_dB)
    sim(ii)=alamouti_sim_mimo(precision,16,Eb_N0_dB(ii));
end

%% ��ͼ
%close(hh);
figure;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2); hold on;
semilogy(Eb_N0_dB,sim,'r*','LineWidth',2); hold on;
axis normal tight; grid on;
legend('L=4 ������','2��2�� �������η���','2��2�� Alamouti����');
xlabel('EbN0 (dB)'); ylabel('BER');
title('16QAM��Beamforming��Alamouti���� (MIMO)');