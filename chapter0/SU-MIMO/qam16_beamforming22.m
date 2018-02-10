% qam16 BER仿真程序
% AWGN的Rayleigh衰落信道MIMO系统仿真（波束成形对比版本）

%% 参数准备
clear all; close all; clc;
N = 2^20;
m = N/4;
maxlen = sqrt(10);
nTx = 2;
nRx = 2;
Eb_N0_dB = -2:1:16;
nErr2 = zeros(1,length(Eb_N0_dB));
theoryBer_nRx22 = zeros(1,length(Eb_N0_dB));
hh = waitbar(1,'请等待，正在计算...');

%% 开始仿真
for ii = 1:length(Eb_N0_dB)
   % 生成qam16信号
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % 使用比特流生成qam16星座图
   s = s./maxlen; % 能量归一化，只改动噪声功率
   n = 1/sqrt(2)*(randn(2,m) + 1j*randn(2,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nTx*nRx,m) + 1j*randn(nTx*nRx,m)); % 瑞利信道
   hEff = h.*exp(-1j*angle(h)); % 发送端波束成形
   sigma=calnoisesigma(16,s,Eb_N0_dB(ii),1); % 噪声功率
   sr = (1/sqrt(nTx))*kron(ones(nTx*nRx,1),s); 
   %生成两根天线的发送序列 成为一个4*N的矩阵
   
   % 模拟经过AWGN瑞利信道
   sray=hEff.*sr;
   y2 = [sray(1,:)+sray(2,:);sray(3,:)+sray(4,:)] + sigma*n; 
   y2 = sum(y2);
   
   % 接收端均衡器
   y2Hat = (maxlen).*y2./sum(hEff,1); 
   %figure; scatter(real(y1Hat),imag(y1Hat));
   
   % 接收端硬判决
   [ip2Hat,uu2] = qam16_demodulation(y2Hat,maxlen);

   % 计算理论线和仿真错误比特数
   nErr2(ii) = size(find(ip- ip2Hat),2);
   theoryBer_nRx22(ii) = berfading(Eb_N0_dB(ii),'qam',16,4);
end
simBer2 = nErr2/N; % 仿真BER

%% Alamouti
precision=2^17;
sim=zeros(1,length(Eb_N0_dB));
for ii=1:length(Eb_N0_dB)
    sim(ii)=alamouti_sim_mimo(precision,16,Eb_N0_dB(ii));
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
title('16QAM的Beamforming和Alamouti仿真 (双发双收MIMO)');