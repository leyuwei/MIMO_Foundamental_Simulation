clc;clear all;close all;

%% 开始仿真Alamouti
ebn0db=-2:1:16;
precision=2^16;
sim=zeros(1,length(ebn0db));
theo=zeros(1,length(ebn0db));
for ii=1:length(ebn0db)
    sim(ii)=alamouti_sim(precision,16,ebn0db(ii));
    theo(ii)=berfading(ebn0db(ii),'qam',16,2);
end

%% 开始仿真Beamforming
maxlen=sqrt(10);
N=precision;
m=N/log2(16);
nTx=2;
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % 使用比特流生成qam16星座图
   s = s./maxlen; % 能量归一化，只改动噪声功率
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % 瑞利信道
   hEff = h.*exp(-1j*angle(h)); % 发送端波束成形
   sigma=calnoisesigma(16,s,ebn0db(ii),1); % 噪声功率
   sr = kron(ones(nTx,1),s)./sqrt(2); %生成两根天线的发送序列 成为一个2*N的矩阵 也可以说是一个二维向量
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   [ip2Hat,uu2] = qam16_demodulation(y2Hat,maxlen);
   nErr2(ii) = size(find(ip- ip2Hat),2);
end
simBer2 = nErr2/N; % 仿真BER （带Beamforming）

%% 绘图
semilogy(ebn0db,theo,'k+-','LineWidth',1);hold on;
semilogy(ebn0db,sim,'r*','LineWidth',2);hold on;
semilogy(ebn0db,simBer2,'b+','LineWidth',2);
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');
title('16QAM的Alamouti仿真');
legend('双径理论值','Alamouti仿真点','Beamforming仿真点');