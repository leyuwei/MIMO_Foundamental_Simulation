% 16QAM
% 2x发射天线 2x单天线用户 MU-MIMO DL BER仿真程序
% 发射端Beamforming，使用ZF-BD算法

%% 参数准备
clear all; close all; clc;
N = 2^17;
m = N / 4;
maxlen = sqrt(10);
nTx = 2;
nRx = 2;
Eb_N0_dB = -2:1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'请等待，正在仿真中...');

%% 开始仿真
for ii = 1:length(Eb_N0_dB)
   % 生成qam信号
   ip = rand(1,N)>0.5;
   ip2 = rand(1,N)>0.5;
   s = qam16_modulation(ip); % 使用比特流调制QAM星座图
   s2 = qam16_modulation(ip2);
   s = s./maxlen; % 能量归一化，只改动噪声功率
   s2 = s2./maxlen;
   n = 1/sqrt(2)*(randn(nRx,1,m) + 1j*randn(nRx,1,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nRx,nTx,m) + 1j*randn(nRx,nTx,m)); % 瑞利信道
   sigma = calnoisesigma(16,s,Eb_N0_dB(ii),1); % 噪声功率
   sr = zeros(nTx,1,m);
   transig = [s;s2];
   for t = 1:nTx
       sr(t,1,:) = (1/sqrt(nTx)) * transig(t,:); 
   end
   %生成两根天线的发送序列 成为一个2*1*m的矩阵
   
   % 模拟经过AWGN瑞利信道
   % CSI分解
   [precode, equalizer, diagonal, sray]=svdChannel(h,sr);
   y = sray + sigma.*n; % sray = H * precode_mat * signal
   
   % 接收端均衡
   y2Hat = zeros(nRx,1,m);
   for jj = 1:m
       y2Hat(:,:,jj) = equalizer(:,:,jj) * y(:,:,jj);
   end
   
   % 接收端信号检测
   outputsig = zeros(nRx,1,m);
   for r = 1:nRx
       outputsig(r,1,:) = y2Hat(r,1,:).*conj(diagonal(r,r,:))./(abs(diagonal(r,r,:))).^2;
   end
   y2Hat = outputsig * sqrt(nTx) * maxlen;
   
   % 接收端硬判决
   ip2Hat1 = qam16_demodulation(y2Hat(1,:),maxlen);
   ip2Hat2 = qam16_demodulation(y2Hat(2,:),maxlen);

   % 计算理论线和仿真错误比特数
   nErr2(ii) = size(find(ip - ip2Hat1),2) + size(find(ip2 - ip2Hat2),2);
   disp(['bit error rate: ' num2str(nErr2(ii)/(nTx*N))]);
end
simBer2 = nErr2 / (2*N); % 仿真BER

%% 绘图
close(hh);
figure;
semilogy(Eb_N0_dB,simBer2,'k*-','LineWidth',1); hold on;
axis normal; grid on;
legend('仿真结果点');
xlabel('EbN0 (dB)'); ylabel('BER');
title('16QAM 2X2 MU-MIMO仿真结果');