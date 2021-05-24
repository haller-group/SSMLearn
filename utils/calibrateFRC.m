function f_red = calibrateFRC(yCal, Omega, V, Tinv, damp, freq)

for iAmp = 1:size(yCal,2)
    zCal = Tinv(V.'*yCal(:,iAmp));
    rhoCal = abs(zCal(1,:));
    f_red(iAmp) = sqrt(rhoCal^2*(freq(rhoCal)-Omega(iAmp))^2+rhoCal^2*(damp(rhoCal))^2);
end

% xCal = {t_sim.', observable(x_sim.')};
% yCal = coordinates_embedding(xCal, SSMDim, 'OverEmbedding', overEmbed);
% etaCal = getProjectedTrajs(yCal, V);
% zCal = transformComplex(Tinv, etaCal);
% rhoCal = mean(abs(zCal{1,2}(1,end-1000:end)));
% fpsi = fsolve(@(fpsi) [damp(rhoCal)*rhoCal + fpsi(1)*sin(fpsi(2)); ...
%     freq(rhoCal) - Omega + fpsi(1)/rhoCal*cos(fpsi(2))],...
%     [0; 0]);
% f_red = abs(fpsi(1));

% calibrationLoads = f_full(end);
% calloadvector = calibrationLoads.*f_vec;
% ICCal = getStaticResponse(K, M, F, calloadvector, 0, 0);
% uCal = observable(ICCal);

% for iAmp = 1:length(amplitudes)
%     [uCal, pos] = max(expAmp{iAmp}(:,2));
%     if iAmp == 4; [uCal, pos] = min(expAmp{iAmp}(:,2)); end
%     Omega = expAmp{iAmp}(pos,1)*7.8;
%     yCal = uCal*V(:,1)./V(floor(embedDim/2)*length(rawColInds)+1,1);
%     zCal = Tinv(V.'*yCal);
%     rhoCal = abs(zCal(1));
%     f_red(iAmp) = sqrt(rhoCal^2*(freq(rhoCal)-Omega)^2+rhoCal^2*(damp(rhoCal))^2);
% end
% f_red = amplitudes*f_red(iAmp)/amplitudes(iAmp)