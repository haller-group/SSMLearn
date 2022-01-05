function BBCInfo = backboneCurves(IMInfo, RDInfo, amplitudeFunction, ...
    maxRho, varargin)
% backboneCurves(IMInfo, RDInfo, amplitudeFunction, maxRho)
% backboneCurves(IMInfo, RDInfo, amplitudeFunction, maxRho, 'Hz')
% backboneCurves(IMInfo, RDInfo, amplitudeFunction, maxRho, 'norm')
% Plot instantaneous amplitude and frequency/damping curves. 
%   
% INPUT
% IMInfo             struct            Manifold.
% RDInfo             struct            Reduced dynamics.
% amplitudeFunction  function handle   Map from y to scalar. Can be
%                                      used to predict forced response 
%                                      for a given component or quantity.
% maxRho             real >0           Normal form maximal plot amplitude.
% scaling            'Hz' or 'norm'    Optionally rescales the frequency
%                                      unit from rad/s to hertz or
%                                      normalized with the eigenfrequency.
% OUTPUT
% BBCInfo            struct            Instantaneous damping, frequency,
%                                      and amplitude arrays parametrized by
%                                      amplitudeNormalForm.
%
% The amplitude is defined as the maximum absolute value reached by the the 
% scalar input function amplitudeFunction along a full rotation in the 
% angle theta for each value of the parametrizing amplitude from 0 to 
% maxRho. If an additional argument is given, the values of frequency and 
% damping are normalized with respect to the linear limit (input 'norm') or 
% in terms of damping ration and Hertz for the frequency.
% For a SSM with more than 1 mode, the backbone curves are computed for 
% each uncoupled limit.

T = RDInfo.transformation.map;
damp = RDInfo.conjugateDynamics.damping;
freq = RDInfo.conjugateDynamics.frequency;
SSM_func = IMInfo.parametrization.map;
k = length(RDInfo.eigenvaluesLinPartFlow);
if nargin(RDInfo.conjugateDynamics.damping) ~= 1 || rem(k,2) ~= 0
    error('Backbone curves not available for this normal form.')
else
    ndof = k/2;
    if length(maxRho) ~= ndof; maxRho = maxRho(1)*ones(1,ndof); end
    
    % Compute instantaneous damping and frequency
    nRhoEvals = 1001;
    rhoMat = zeros(ndof,nRhoEvals); instDampCurve = rhoMat;
    instFreqCurve = rhoMat; instAmp = rhoMat;
    for idof = 1:ndof
        rhoMat(idof,:) = linspace(0,maxRho(idof),nRhoEvals);
        rhoPlot = zeros(ndof,nRhoEvals); rhoPlot(idof,:) = rhoMat(idof,:);
        dampTemp = damp(rhoPlot); freqTemp = freq(rhoPlot);
        instDampCurve(idof,:) = dampTemp(idof,:);
        instFreqCurve(idof,:) = freqTemp(idof,:);
        % Compute instantaneous amplitude
        theta_plot = linspace(0,2*pi,51); theta_plot = theta_plot(1:end-1);
        [RR,TT] = meshgrid(rhoPlot(idof,:),theta_plot);
        ZZ = RR.*exp(1i*TT);
        zEval = zeros(ndof,numel(RR)); zEval(idof,:) = transpose(ZZ(:));
        zEval = transformationComplexConj(zEval);
        ampEval = abs(amplitudeFunction( SSM_func( T(zEval) ) ));
        instAmp(idof,:) = max(reshape(transpose(ampEval),size(RR)),[],1);
    end
    
    % Amplitude as maximum value along theta
    BBCInfo = struct('damping',instDampCurve,'frequency',instFreqCurve,...
        'amplitude',instAmp,'amplitudeNormalForm',rhoMat);
    for idof = 1:ndof
        % Plot results
        fig = customFigure('subPlot',[1 2]);
        if isempty(varargin)==0
            if strcmp(varargin{:},'Hz')==1
                subplot(121); hold on; grid on; box on;
                plot(-instDampCurve(idof,:)./instFreqCurve(idof,:)*100, ...
                    instAmp(idof,:),'Linewidth',2)
                if ndof == 1
                    xlabel('Damping ratio [%]')
                else
                    xlabel(['$c_{' num2str(idof) '}/\omega_{' ...
                       num2str(idof) '} \, [$\%$]$'],'Interpreter','latex')
                end
                subplot(122); hold on; grid on; box on;
                plot(instFreqCurve(idof,:)/2/pi,instAmp(idof,:),...
                                                             'Linewidth',2)
                if ndof == 1
                    xlabel('Frequency [Hz]')
                else
                    xlabel(['$\omega_{' num2str(idof) ...
                             '}/(2\pi) \, [$Hz$]$'],'Interpreter','latex')
                end
            else
                subplot(121); hold on; grid on; box on;
                plot(instDampCurve(idof,:)./instDampCurve(idof,1),...
                    instAmp(idof,:),'Linewidth',2)
                if ndof == 1
                    xlabel('$c/c(0)$','Interpreter','latex')
                else
                    xlabel(['$c_{' num2str(idof) '}/c_{' num2str(idof) ...
                        '}(0)$'],'Interpreter','latex')
                end
                subplot(122); hold on; grid on; box on;
                plot(instFreqCurve(idof,:)/instFreqCurve(idof,1),...
                    instAmp(idof,:),'Linewidth',2)
                if ndof == 1
                    xlabel('$\omega/\omega(0)$','Interpreter','latex')
                else
                    xlabel(['$\omega_{' num2str(idof) '}/\omega_{' ...
                        num2str(idof) '}(0)$'],'Interpreter','latex')
                end
            end
        else
            subplot(121); hold on; grid on; box on;
            plot(instDampCurve(idof,:),instAmp(idof,:),'Linewidth',2)
            if ndof == 1
                xlabel('Damping [1/s]')
            else
                xlabel(['$c_{' num2str(idof) '}$'],'Interpreter','latex')
            end
            subplot(122); hold on; grid on; box on;
            plot(instFreqCurve(idof,:),instAmp(idof,:),'Linewidth',2)
            xlabel('Frequency [rad/s]')
            if ndof == 1
                xlabel('Frequency [rad/s]')
            else
                xlabel(['$\omega_{' num2str(idof) ...
                                  '} \, [$rad/s$]$'],'Interpreter','latex')
            end
        end
        subplot(121);
        ylabel('Amplitude')
        subplot(122);
        ylabel('Amplitude')
    end
end
end