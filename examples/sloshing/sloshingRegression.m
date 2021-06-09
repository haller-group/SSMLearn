%% Finding a 2D SSM from sloshing data
% 

clearvars
close all
clc

%% Example setup
frcdir = 'frcdata/';
expAmp{1} = load([frcdir,'Measurements_A=0.09%.txt']);
expAmp{2} = load([frcdir,'Measurements_A=0.17%.txt']);
expAmp{3} = load([frcdir,'Measurements_A=0.32%.txt']);
% expAmp{4} = load([frcdir,'Measurements_A=0.64%.txt']);
amplitudes = [0.09 0.17, 0.32];
% amplitudes = [0.09 0.17, 0.32, 0.64];

%% Regression
Exp = vertcat(expAmp{1,:});
indexAmplitudes = [];
for i = 1:length(amplitudes)
    indexAmplitudes = [indexAmplitudes; i*ones(size(expAmp{i},1),1)];
end
amplitudes = [0.09 0.17, 0.32, 0.64]
ExpOmega = Exp(:,1)*7.8;
ExpX = Exp(:,2);
% opts = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1e4, 'MaxFunctionEvaluations', 1e5, 'SpecifyObjectiveGradient', true);
% func = @(z) funct(z, ExpOmega, ExpX, indexAmplitudes)
% zOpt = fsolve(func, [-0.0545, 7.79, 0, -0.0158, 0.1, 0.2, 0.36], opts)
opts = optimoptions('fminunc', 'Display', 'iter', 'MaxIterations', 1e4, 'MaxFunctionEvaluations', 1e5, 'SpecifyObjectiveGradient', true);
func = @(z) functmin(z, ExpOmega, ExpX, indexAmplitudes);
zOpt = fminunc(func, [-0.0543, 7.7930, -0.0005, -0.0158, 0.1043, 0.1929, 0.3621, 0.58], opts)
% zOpt = [-0.0543, 7.7930, -0.0005, -0.0158, 0.1043, 0.1929, 0.3621, 0.64]
damp = @(x) zOpt(1) + zOpt(3)*x.^2;
freq = @(x) zOpt(2) + zOpt(4)*x.^2;
f_red = zOpt(4+1:4+3);
f_red(4) = 0.64
expAmp{4} = load([frcdir,'Measurements_A=0.64%.txt']);

%%
w_span = [0.77, 1.06]*7.8;

% Compute FRC analytically
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
for iAmp = 1:length(amplitudes)
%     rhoTip = max(ExpX);
    rhoTip = abs(fsolve(@(rho) 1e5*(f_red(iAmp)-(rho*damp(rho))), max(ExpX), options));
    rhoFRC = logspace(log10(rhoTip*0.03), log10(rhoTip), 300);
    rhoFRC = [rhoFRC, -fliplr(rhoFRC)];
    OmegaFRC(iAmp,:) = 1/7.8*real(freq(rhoFRC) + -1./rhoFRC.*sqrt(f_red(iAmp)^2-(rhoFRC.*damp(rhoFRC)).^2));
    uFRC(iAmp,:) = abs(rhoFRC);
end


%% Plot
figure; hold on; colors = colororder;
% plot(frq, amp,'k','DisplayName', 'Backbone - SSMlearn')
for iAmp = 1:length(amplitudes)
    plot(OmegaFRC(iAmp,:), uFRC(iAmp,:), 'LineWidth', 2, 'Color', colors(iAmp+1,:),...
        'DisplayName', ['SSMLearn A = ', num2str(amplitudes(iAmp)), ' %'])
    plot(expAmp{iAmp}(:,1), expAmp{iAmp}(:,2), '.', 'MarkerSize', 12, ...
        'Color', colors(iAmp+1,:), 'DisplayName', ['Exp. A = ', num2str(amplitudes(iAmp)), ' %'])
end
xlim(w_span/7.8)
ylim([0,10])
xlabel('Excitation frequency $\Omega$ [normalized]', 'Interpreter', 'latex')
ylabel('Amplitude $\hat{X}$ (\%)', 'Interpreter', 'latex')
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)
legend

function [f, Df] = funct(z, O, X, indexAmplitudes)
ord = 4;
a = z(1)*X + z(3)*X.^3;
b = z(2) + z(4)*X.^2;
f = (O - b).^2.*X.^2 - (z(ord+indexAmplitudes)').^2 + a.^2;
Dfab = [2*a.*X, -2*(O - b).*X.^2, 2*a.*X.^3, -2*(O - b).*X.^4];
Dff = [-2*z(ord+1).*(indexAmplitudes==1), -2*z(ord+2).*(indexAmplitudes==2), -2*z(ord+3).*(indexAmplitudes==3)];%, -2*z(ord+4).*(indexAmplitudes==4)];
Df = [Dfab, Dff];
end

function [f, Df] = functmin(z, O, X, indexAmplitudes)
ord = 4;
a = z(1)*X + z(3)*X.^3;
b = z(2) + z(4)*X.^2;
P = (O - b).^2.*X.^2 - (z(ord+indexAmplitudes)').^2 + a.^2;
f = P'*P;
Dfab = [2*a.*X, -2*(O - b).*X.^2, 2*a.*X.^3, -2*(O - b).*X.^4];
Dff = [-2*z(ord+1).*(indexAmplitudes==1), -2*z(ord+2).*(indexAmplitudes==2), -2*z(ord+3).*(indexAmplitudes==3), -2*z(ord+4).*(indexAmplitudes==4)];
Df = 2*P'*[Dfab, Dff];
end
