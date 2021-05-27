clear all
close all
clc

a0 = -0.061747;
a2 = -0.070295;
w0 = 7.8098;
w2 = -1.8774;
ts = 0.01;
damp_an = @(r) a0+a2*r.^2;
freq_an = @(r) (w0+w2*r.^2)/2/pi;

f = @(t,x) [a0*x(1) - w0*x(2) + (x(1).^2+x(2).^2).*(a2*x(1) - w2*x(2)); ...
            a0*x(2) + w0*x(1) + (x(1).^2+x(2).^2).*(a2*x(2) + w2*x(1))];

[t,x] = ode45(f,[0:ts:80],[1; 0]);
x = x(:,1);
figure(1); clf; hold on; grid on; box on;
plot(t,x,'Linewidth',2)
xlabel('time')
ylabel('signal')
set(gca,'fontname','helvetica')
set(gca,'fontsize',16)
%%
[amp,freq,damp,time] = PFF(t,x);

plot(time,amp,'r','Linewidth',2)
plot(time,-amp,'r','Linewidth',2)

figure(2); clf;
subplot(121); hold on; grid on; box on;
r_disp = linspace(0,1,51);
plot(freq_an(r_disp),r_disp,'Linewidth',2)
plot(freq,amp,'Linewidth',2)
xlabel('frequency [Hz]')
ylabel('amplitude')
set(gca,'fontname','helvetica')
set(gca,'fontsize',16)
subplot(122); hold on; grid on; box on;
r_disp = linspace(0,1,51);
plot(damp_an(r_disp),r_disp,'Linewidth',2)
plot(damp,amp,'Linewidth',2)
xlabel('damping [1/s]')
ylabel('amplitude')
set(gca,'fontname','helvetica')
set(gca,'fontsize',16)




