% QPSK BER仿真程序
% AWGN的Rayleigh衰落信道MISO系统仿真（波束成形对比版本）

%% 参数准备
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
hh = waitbar(1,'请等待，正在计算...');

%% 开始仿真
for ii = 1:length(Eb_N0_dB)
   % 生成qpsk信号
   ip = rand(1,N)>0.5;
   s = qpsk_modulation(ip); % 使用比特流生成qpsk星座图
   s = s./maxlen; % 能量归一化，只改动噪声功率
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % 瑞利信道
   hEff = h.*exp(-1j*angle(h)); % 发送端波束成形
   sigma=calNoisePower(4,s,Eb_N0_dB(ii)); % 噪声功率
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %生成两根天线的发送序列 成为一个2*N的矩阵 也可以说是一个二维向量
   
   % 模拟经过AWGN瑞利信道
   y1 = sum(h.*sr,1) + sigma*n; 
   y2 = sum(hEff.*sr,1) + sigma*n; 

   % 接收端均衡器
   y1Hat = maxlen.*y1./sum(h,1); 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   % figure; scatter(real(y1Hat),imag(y1Hat));
   
   % 接收端硬判决
   [ip1Hat,uu] = qpsk_demodulation(y1Hat);
   [ip2Hat,uu2] = qpsk_demodulation(y2Hat);

   % 计算理论线和仿真错误比特数
   nErr1(ii) = size(find(ip- ip1Hat),2);
   nErr2(ii) = size(find(ip- ip2Hat),2);
   theoryBer_nRx2(ii) = berfading(Eb_N0_dB(ii),'psk',4,1);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'psk',4,2);
end
simBer1 = nErr1/N; % 仿真BER （不带Beamforming）
simBer2 = nErr2/N; % 仿真BER （带Beamforming）

%% 绘图
close(hh);
figure;
semilogy(Eb_N0_dB,theoryBer_nRx2,'b*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(Eb_N0_dB,simBer1,'r*','LineWidth',2); hold on;
semilogy(Eb_N0_dB,simBer2,'mx','LineWidth',2);
axis normal tight; grid on;
legend('L=1 理论线','L=2 理论线','2发1收 无波束成形仿真','2发1收 波束成形仿真');
xlabel('EbN0 (dB)'); ylabel('BER');
title('QPSK在AWGN瑞利衰落信道下的MISO系统仿真 (波束成形)');