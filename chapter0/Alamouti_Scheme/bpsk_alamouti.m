% BPSK BER仿真程序
% AWGN的Rayleigh衰落信道MISO系统仿真（Alamouti对比版本）

%% 参数准备
clear all; close all; clc; 
N = 2^16;
m = N;
nTx = 2;
maxlen = 1;
ebn0db = -2:16;
Xs = [-1 1]; % 星座点，供ML判决使用
BER = zeros(1,length(ebn0db));
theoryBer_nRx22 = zeros(1,length(ebn0db));
% 生成bpsk序列
ip = rand(1,N)>0.5;
ST = 2*ip-1;
hh = waitbar(0,'请稍等 正在计算...');

%% 开始仿真Alamouti
for n=1:length(ebn0db)
    waitbar(n/length(ebn0db),hh,['Alamouti计算：' num2str(100*n/length(ebn0db)) '%']);
    ST_estimated=zeros(1,length(ST));
    for k=1:2:length(ST)
        % 对原始信号进行串并变换，分流成两路
        s0=ST(k);
        s1=ST(k+1);
        
        % 初始化噪声
        sigma=calNoisePower(2,ST,ebn0db(n));
        Noise0=sigma.*(randn+1j*randn)./sqrt(2);
        Noise1=sigma.*(randn+1j*randn)./sqrt(2);
        
        % 初始化两径iid的瑞利衰落信道系数
        h0=(randn+1j*randn)*sqrt(1/2);
        h1=(randn+1j*randn)*sqrt(1/2);
        hEq=(abs(h0))^2+(abs(h1))^2;

        % 根据Alamouti原始论文 式(11) 接收信号表示如下
        r0=h0*s0+h1*s1+Noise0;  
        r1=-h0*conj(s1)+h1*conj(s0)+Noise1;      

        % 根据Alamouti原始论文 式(12) 合并模式表示如下
        S0=conj(h0).*r0+h1.*conj(r1);
        S1=conj(h1).*r0-h0.*conj(r1);
        
        % 基于ML准则的判决
        S=[S0,S1]./hEq;
        % 根据ML准则的要求，计算Squared Euclidean Distance即平方欧氏距离
        % 简单来说就是判别其距离哪一个坐标点更近
        d11=(real(S(1))+1)^2+(imag(S(1))-0)^2;
        d12=(real(S(1))-1)^2+(imag(S(1))-0)^2;
        d21=(real(S(2))+1)^2+(imag(S(2))-0)^2;
        d22=(real(S(2))-1)^2+(imag(S(2))-0)^2;
        % 根据ML准则的判决要求，找出最小距离
        [ds1_min,position1]=min([d11,d12]);
        [ds2_min,position2]=min([d21,d22]);
        % 解码
        s1_estemated=Xs(position1);         
        s2_estemated=Xs(position2);
        ST_estimated(k)=s1_estemated;
        ST_estimated(k+1)=s2_estemated; 
    end
    % 解调
    R=sign(real(ST_estimated));
    % 计算BER
    no_errors=size(find(ST-R),2);
    BER(1,n)=no_errors/length(ST);
    theoryBer_nRx22(n) = berfading(ebn0db(n),'psk',2,2);
end

%% 开始仿真Beamforming
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = 2*ip-1;
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % 瑞利信道
   sigma=calNoisePower(2,s,ebn0db(ii)); % 噪声功率
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %生成两根天线的发送序列 成为一个2*N的矩阵 也可以说是一个二维向量
   hEff = h.*exp(-1j*angle(h));
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = y2./sum(hEff,1); 
   ip2Hat = real(y2Hat)>0;
   nErr2(ii) = size(find([ip- ip2Hat]),2);
end
simBer2 = nErr2/N; % 仿真BER （带Beamforming）

%% 绘图
close(hh);figure;
semilogy(ebn0db,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(ebn0db,simBer2,'bx','LineWidth',2); hold on;
semilogy(ebn0db,BER,'r*','LineWidth',2);
legend('L=2 理论线','2发1收 Beamforming仿真','2发1收 Alamouti仿真');
title('BPSK在AWGN的瑞利衰落信道下的MISO系统仿真（Alamouti）');
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');