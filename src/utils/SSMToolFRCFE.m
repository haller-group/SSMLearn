function [FRCSSMTool] = SSMToolFRCFE(M,C,K,fnl,fext,outdof,epsilon,masterModes,order,freqRange,mFreqs,oid)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Create
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Forcing assumped single harmonic cosine
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas,epsilon(1));

% SSM setup
S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex');
set(S.FRCOptions, 'nCycle',500, 'initialSolver', 'forward');
set(S.contOptions, 'PtMX', 300, 'h_max', 0.5);
set(S.FRCOptions, 'omegaSampStyle', 'cocoBD');

% FRC computation
FRCSSMTool = struct();
if isscalar(epsilon)
    FRC = S.SSM_isol2ep(oid,masterModes, order, mFreqs, 'freq', freqRange,outdof);
    ampFRCoutTemp = [FRC.Aout_frc];
    if size(ampFRCoutTemp,2) > 1
        ampFRCout = zeros(1,size(ampFRCoutTemp,1),size(ampFRCoutTemp,2));
        ampFRCout(1,:,:) = ampFRCoutTemp;
    else
        ampFRCout = transpose(ampFRCoutTemp);
    end
    FRCSSMTool.(['F1']) = struct('Freq',transpose([FRC.om]),'Amp',...
        ampFRCout,'Nf_Amp',transpose([FRC.rho]),'Nf_Phs',transpose([FRC.th]),'Stab',transpose([FRC.st]));
else
    FRC = S.SSM_epSweeps(oid,masterModes,order,mFreqs,epsilon,freqRange,outdof);
    for ii = 1:length(epsilon)
        if isempty(FRC.FRCom)
            iFRC = struct('om',[],'Aout_frc',[],'rho',[],'th',[],'st',[]);
            ampFRCout = [];
        else
        iFRC = FRC.FRCom{ii};
        ampFRCoutTemp = [iFRC.Aout_frc];
        if size(ampFRCoutTemp,2) > 1
            ampFRCout = zeros(1,size(ampFRCoutTemp,1),size(ampFRCoutTemp,2));
            ampFRCout(1,:,:) = ampFRCoutTemp;
        else
            ampFRCout = transpose(ampFRCoutTemp);
        end
        end
        FRCSSMTool.(['F' num2str(ii)]) = struct('Freq',transpose([iFRC.om]),'Amp',...
            ampFRCout,'Nf_Amp',transpose([iFRC.rho]),'Nf_Phs',transpose([iFRC.th]),'Stab',transpose([iFRC.st]));
    end
end
end

