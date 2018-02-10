% qam16 BER仿真程序
% AWGN的Rayleigh衰落信道MISO系统仿真（Alamouti对比版本）

%% 参数准备
clear all; close all; clc; 
N = 2^16;
m = N/4;
nTx = 2;
maxlen = sqrt(10);
ebn0db = -2:16;
Xs = [1+1j 3+1j 3+3j 1+3j -1+3j -1+1j -3+3j -3+1j -3-1j -3-3j -1-3j -1-1j 1-1j 3-1j 1-3j 3-3j]; % 星座点，供ML判决使用
BER = zeros(1,length(ebn0db));
theoryBer_nRx22 = zeros(1,length(ebn0db));
% 生成qam16序列
ip = rand(1,N)>0.5;
ST = qam16_modulation(ip);
ST=ST./maxlen;
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
        sigma=calNoisePower(16,ST,ebn0db(n));
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
        Xss=[real(Xs);imag(Xs)];
        d11=squaredis(S(1),[Xss(1,1) Xss(2,1)],maxlen);
        d12=squaredis(S(1),[Xss(1,2) Xss(2,2)],maxlen);
        d13=squaredis(S(1),[Xss(1,3) Xss(2,3)],maxlen);
        d14=squaredis(S(1),[Xss(1,4) Xss(2,4)],maxlen);
        d15=squaredis(S(1),[Xss(1,5) Xss(2,5)],maxlen);
        d16=squaredis(S(1),[Xss(1,6) Xss(2,6)],maxlen);
        d17=squaredis(S(1),[Xss(1,7) Xss(2,7)],maxlen);
        d18=squaredis(S(1),[Xss(1,8) Xss(2,8)],maxlen);
        d19=squaredis(S(1),[Xss(1,9) Xss(2,9)],maxlen);
        d110=squaredis(S(1),[Xss(1,10) Xss(2,10)],maxlen);
        d111=squaredis(S(1),[Xss(1,11) Xss(2,11)],maxlen);
        d112=squaredis(S(1),[Xss(1,12) Xss(2,12)],maxlen);
        d113=squaredis(S(1),[Xss(1,13) Xss(2,13)],maxlen);
        d114=squaredis(S(1),[Xss(1,14) Xss(2,14)],maxlen);
        d115=squaredis(S(1),[Xss(1,15) Xss(2,15)],maxlen);
        d116=squaredis(S(1),[Xss(1,16) Xss(2,16)],maxlen);
        d21=squaredis(S(2),[Xss(1,1) Xss(2,1)],maxlen);
        d22=squaredis(S(2),[Xss(1,2) Xss(2,2)],maxlen);
        d23=squaredis(S(2),[Xss(1,3) Xss(2,3)],maxlen);
        d24=squaredis(S(2),[Xss(1,4) Xss(2,4)],maxlen);
        d25=squaredis(S(2),[Xss(1,5) Xss(2,5)],maxlen);
        d26=squaredis(S(2),[Xss(1,6) Xss(2,6)],maxlen);
        d27=squaredis(S(2),[Xss(1,7) Xss(2,7)],maxlen);
        d28=squaredis(S(2),[Xss(1,8) Xss(2,8)],maxlen);
        d29=squaredis(S(2),[Xss(1,9) Xss(2,9)],maxlen);
        d210=squaredis(S(2),[Xss(1,10) Xss(2,10)],maxlen);
        d211=squaredis(S(2),[Xss(1,11) Xss(2,11)],maxlen);
        d212=squaredis(S(2),[Xss(1,12) Xss(2,12)],maxlen);
        d213=squaredis(S(2),[Xss(1,13) Xss(2,13)],maxlen);
        d214=squaredis(S(2),[Xss(1,14) Xss(2,14)],maxlen);
        d215=squaredis(S(2),[Xss(1,15) Xss(2,15)],maxlen);
        d216=squaredis(S(2),[Xss(1,16) Xss(2,16)],maxlen);
        % 根据ML准则的判决要求，找出最小距离
        [ds1_min,position1]=min([d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111,d112,d113,d114,d115,d116]);
        [ds2_min,position2]=min([d21,d22,d23,d24,d25,d26,d27,d28,d29,d210,d211,d212,d213,d214,d215,d216]);
        % 解码
        s1_estemated=Xs(position1);         
        s2_estemated=Xs(position2);
        ST_estimated(k)=s1_estemated;
        ST_estimated(k+1)=s2_estemated; 
    end
    % 解调
    [R,~]=qam16_demodulation(ST_estimated,maxlen);
    % 计算BER
    no_errors=size(find(ip-R),2);
    BER(1,n)=no_errors/length(ip);
    theoryBer_nRx22(n) = berfading(ebn0db(n),'qam',16,2);
end

%% 开始仿真Beamforming
nErr2 = zeros(1,length(ebn0db));
for ii = 1:length(ebn0db)
   ip = rand(1,N)>0.5;
   s = qam16_modulation(ip); % 使用比特流生成qam16星座图
   s = s./maxlen; % 能量归一化，只改动噪声功率
   n = 1/sqrt(2)*(randn(1,m) + 1j*randn(1,m)); % 高斯白噪声
   h = 1/sqrt(2)*(randn(nTx,m) + 1j*randn(nTx,m)); % 瑞利信道
   hEff = h.*exp(-1j*angle(h)); % 发送端波束成形
   sigma=calNoisePower(16,s,ebn0db(ii)); % 噪声功率
   sr = (1/sqrt(nTx))*kron(ones(nTx,1),s); %生成两根天线的发送序列 成为一个2*N的矩阵 也可以说是一个二维向量
   y2 = sum(hEff.*sr,1) + sigma*n; 
   y2Hat = maxlen.*y2./sum(hEff,1); 
   [ip2Hat,uu2] = qam16_demodulation(y2Hat,maxlen);
   nErr2(ii) = size(find(ip- ip2Hat),2);
end
simBer2 = nErr2/N; % 仿真BER （带Beamforming）

%% 绘图
close(hh);figure;
semilogy(ebn0db,theoryBer_nRx22,'k*-','LineWidth',1); hold on;
semilogy(ebn0db,simBer2,'bx','LineWidth',2); hold on;
semilogy(ebn0db,BER,'r*','LineWidth',2);
legend('L=2 理论线','2发1收 Beamforming仿真','2发1收 Alamouti仿真');
title('16QAM在AWGN的瑞利衰落信道下的MISO系统仿真（Alamouti）');
axis normal tight; grid on;
xlabel('Eb/N0 (dB)'); ylabel('BER');