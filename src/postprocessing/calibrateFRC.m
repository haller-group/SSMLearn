function fRed = calibrateFRC(IMInfo, RDInfo, yCal, Omega)
%   fRed = calibrateFRC(IMInfo, RDInfo, yCal, Omega)
%   Compute the forcing amplitudes in the normal form such that the FRC
%   gives the responses yCal at forcing frequencies Omega
%   
%   INPUT
%   IMInfo  struct
%   RDInfo  struct
%   yCal    (n x nAmp)     calibration amplitudes in observable space
%   Omega   (nAmp x 1)     calibration frequencies

SSMChart = IMInfo.chart.map;
damp = RDInfo.conjugateDynamics.damping;
freq = RDInfo.conjugateDynamics.frequency;
invT = RDInfo.inverseTransformation.map;
ndof = IMInfo.parametrization.dimension/2;

nAmp = size(yCal,2);
fRed = zeros(nAmp,ndof);
for iAmp = 1:nAmp
    zCal = invT(SSMChart(yCal(:,iAmp)));
    rhoCal = abs(zCal(1:ndof,:));
    fRed(iAmp,:) = sqrt(rhoCal.^2.*(freq(rhoCal)-Omega(iAmp)).^2 + rhoCal.^2.*(damp(rhoCal)).^2);
end
