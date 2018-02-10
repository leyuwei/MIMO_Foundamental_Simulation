% Alamouti MIMO ����

clc;clear all;close all;
ebn0db=-2:1:16;

%% 16QAM
precision=2^16;
sim=zeros(1,length(ebn0db));
theo=zeros(1,length(ebn0db));
h=waitbar(1,'����16QAM');
for ii=1:length(ebn0db)
    sim(ii)=alamouti_sim_mimo(precision,16,ebn0db(ii));
    theo(ii)=berfading(ebn0db(ii),'qam',16,4);
end
close(h);
semilogy(ebn0db,theo,'k+-','LineWidth',1);hold on;
semilogy(ebn0db,sim,'r*','LineWidth',2);
axis normal tight; grid on;
title('16QAM��˫��˫��Alamouti����');
legend('L=4 ����ֵ','Alamouti�����');

%% QPSK
figure;
precision=2^19;
sim=zeros(1,length(ebn0db));
theo=zeros(1,length(ebn0db));
h=waitbar(1,'����QPSK');
for ii=1:length(ebn0db)
    sim(ii)=alamouti_sim_mimo(precision,4,ebn0db(ii));
    theo(ii)=berfading(ebn0db(ii),'psk',4,4);
end
close(h);
semilogy(ebn0db,theo,'k+-','LineWidth',1);hold on;
semilogy(ebn0db,sim,'r*','LineWidth',2);
axis normal tight; grid on;
title('QPSK��˫��˫��Alamouti����');
legend('L=4 ����ֵ','Alamouti�����');

%% BPSK
figure;
precision=2^21;
sim=zeros(1,length(ebn0db));
theo=zeros(1,length(ebn0db));
h=waitbar(1,'����BPSK');
for ii=1:length(ebn0db)
    sim(ii)=alamouti_sim_mimo(precision,2,ebn0db(ii));
    theo(ii)=berfading(ebn0db(ii),'psk',2,4);
end
close(h);
semilogy(ebn0db,theo,'k+-','LineWidth',1);hold on;
semilogy(ebn0db,sim,'r*','LineWidth',2);
axis normal tight; grid on;
title('BPSK��˫��˫��Alamouti����');
legend('L=4 ����ֵ','Alamouti�����');