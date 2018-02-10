% bpsk BER仿真程序
% AWGN的Rayleigh衰落信道MIMO系统仿真（波束成形对比版本）

%% 参数准备
clear all; close all; clc;
N = 2^16;
m = N;
maxlen = 1;
nTx = 2;
nRx = 2;
Eb_N0_dB = 1000:-1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'请等待，正在计算...');

%% 开始仿真
for ii = 1:length(Eb_N0_dB)
   % 生成bpsk信号
   ip = rand(1,N)>0.5;
   s = ip*2-1; % 使用比特流生成bpsk星座图
   s = s./maxlen; % 能量归一化，只改动噪声功率
   n = 1/sqrt(2)*(randn(nRx,nRx,m) + 1j*randn(nRx,nRx,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % 瑞利信道
   sigma=calnoisesigma(2,s,Eb_N0_dB(ii),1); % 噪声功率
   sr = zeros(nTx,nRx,m);
   for r=1:nRx
       for t=1:nTx
           sr(t,r,:)=(1/sqrt(nTx))*s; 
       end
   end
   %生成两根天线的发送序列 成为一个2*2*m的矩阵
   
   % 模拟经过AWGN瑞利信道
   % CSI分解
   [precode, equalizer, diagonal, sray]=svdChannel(h,sr);
   y = sray + sigma.*n; 
   
   % 接收端均衡器
   y2Hat=zeros(nRx,nRx,m);
   for jj=1:m
       y2Hat(:,:,jj) = equalizer(:,:,jj)*y(:,:,jj)./abs(sum(diagonal(:,:,jj)));
   end
   
   % 接收端MIMO信号检测
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
   recvsig=tz_rr(1,:)+tz_rr(2,:);
   y2Hat2 = recvsig./sqrt(nTx);
   y2Hat2 = reshape(y2Hat2,1,m);
   
   % 接收端硬判决
   ip2Hat = y2Hat2>0;

   % 计算理论线和仿真错误比特数
   nErr2(ii) = size(find(ip- ip2Hat),2);
   disp(['error rate: ' num2str(nErr2(ii)/N)]);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',2,4);
end
simBer2 = nErr2/N; % 仿真BER

%% Alamouti
precision=2^21;
sim=zeros(1,length(Eb_N0_dB));
for ii=1:length(Eb_N0_dB)
    sim(ii)=alamouti_sim_mimo(precision,2,Eb_N0_dB(ii));
end

%% 绘图
close(hh);
figure;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2); hold on;
semilogy(Eb_N0_dB,sim,'r*','LineWidth',2); hold on;
axis normal tight; grid on;
legend('L=4 理论线','2发2收 波束成形仿真','2发2收 Alamouti仿真');
xlabel('EbN0 (dB)'); ylabel('BER');
title('BPSK的Beamforming和Alamouti仿真 (MIMO)');