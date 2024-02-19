```
matlab

%%
% 数字基带传输和接收
clear all;
close all;

N = 100;
N_sample = 8; % 每码元抽样点数
Ts = 1;
dt = Ts/N_sample;
t = 0:dt:(N*N_sample-1)*dt;

gt = ones(1,N_sample); % 数字基带波形
d = sign(randn(1,N)); % 输入数字序列
a = sigexpand(d,N_sample);
st = conv(st,ht1);

ht1 = gt;
rt1 = conv(st,ht1);
ht2 = 5*sinc(5*(t-5)/Ts);
rt2 = conv(st,ht2);

figure(1)
subplot(321)
plot(t,st(1:length(t)));
axis([0 20 -1.5 1.5]);
ylabel('输入双极性NRZ数字基带波形');

subplot(322)
stem(t,a);
axis([0 20 -1.5 1.5]);
ylabel('输入数字序列');

subplot(323)
plot(t,[0 rt1(1:length(t)-1)]/8);
axis([0 20 -1.5 1.5]);
ylabel('方波滤波后的输出');


%不想写啦，数字基带用处本来就不大，我要去写数字频带了

```
