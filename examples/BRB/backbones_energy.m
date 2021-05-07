% Need to simulate with different sampling time
dtAcc = aData{1,1}(2)-aData{1,1}(1);
dtDIC = uData{1,1}(2)-uData{1,1}(1);
if dtAcc~=dtDIC
    uDataRef=cell(1,2);
    tDis = uData{1,1}(1):dtAcc:uData{1,1}(end); tDis = tDis-tDis(1);
    
    [tDis,zRecRef] = ode45(@(t,x) N(x),tDis,zData{1,2}(:,1));
    yRecRef = SSMFunction(T(transpose(zRecRef)));
    uDataRef{1,1} = transpose(tDis); uDataRef{1,2} = yRecRef(1:206,:);
end
rhoDis = transpose(abs(zRecRef(:,1)));

% Numerical Differentiation on the reconstructed data
tRec = uDataRef{1,1};
URec = uDataRef{1,2}/1e3; n_steps = 3;
[VRec,URec,tRec] = finitetimedifference(URec,tRec,n_steps);
rhoRec = rhoDis(1*n_steps+1:end-1*n_steps);
energyDIC = 0.5*mean(VRec.^2,1)*1.796;
[energyDIC,locs] = findpeaks(energyDIC);
rhoRec = rhoRec(locs);
dampDIC = damp(rhoRec); freqDIC = freq(rhoRec);
damprDIC = -dampDIC./freqDIC; freqDIC = freqDIC/2/pi;


figure; clf;
subplot(121); hold on; box on; grid on;
plot(damprDIC*100,energyDIC*1000,'Linewidth',2)
xlabel('$\xi \, [$\%$]$','interpreter','latex'); 
ylabel('$K \, [$mJ$]$','interpreter','latex'); 
set(gca,'fontname', fontName); set(gca,'fontsize', fontSize); 
xlim([0.1 0.6])
ylim([0.02 20])
set(gca,'yscale','log')
subplot(122); hold on; box on; grid on;
plot(freqDIC,energyDIC*1000,'Linewidth',2)
xlabel('$\omega \, [$Hz$]$','interpreter','latex'); 
ylabel('$K \, [$mJ$]$','interpreter','latex'); 
set(gca,'fontname', fontName); set(gca,'fontsize', fontSize); 
xlim([80 80.4])
ylim([0.02 20])
set(gca,'yscale','log')
