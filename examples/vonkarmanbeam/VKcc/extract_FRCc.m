function FRC = extract_FRCc(obj, parName, parRange, ORDER)
%  EXTRACT_FRC This function extracts the forced response curve (FRC) for
%  systems that may have internal resonances. The FRC computation is based
%  on SSM computation. An appropriate SSM is constructed based on the
%  resonant spectrum of system. The FRC is computed for the reduced system
%  either by level-set technique or continuation method. The
%  obtained FRC is finally mapped back to physical coordinates. By default,
%  level-set method is used for 2D, underdamped SSMs. In any case,
%  the continuation method is employed for higher dimensional SSMs.
%
% FRC = EXTRACT_FRC(OBJ, PARNAME, PARRANGE, ORDER)
%
% parName:  'amp'/'freq' continuation parameter. If parName='amp', FRC is
%           obtained with respect to forcing amplitude epsilon.
%           Specifically, SSM is constructed with assumption that Omega is
%           equal to obj.System.Omega. It follows that the constructed SSM
%           could have dimension higher than two. Then FRC is obtained by
%           level set method or parameter continuation in epsilon.
%           If parName='freq', FRC is obtained with respect to forcing
%           excitation frequency. In this case, SSM is constructed based on
%           the info of parRange. The range of frequency is divided into
%           some subintervals with a single natural frequency is included
%           in each subinterval, and the center point of two adjacent
%           natural frequencies defines the boundary point of two
%           consecutive subintervals. Within each subinterval, appropriate
%           SSM is constructed and then FRC is obtained.
% parRange: continuation domain of parameter
% order:    order of SSM expansion to be used for FRC computation
f1 = figure(1); hold on;
f2 = figure(2); hold on;
figs = [f1, f2];
colors = get(0,'defaultaxescolororder');
totalComputationTime = zeros(size(ORDER));

for j = 1:numel(ORDER)
    order = ORDER(j);
    startFRC = tic;
    if isempty(obj.System.spectrum)
        [~,~,~] = obj.System.linear_spectral_analysis();
    end
    lambda  = obj.System.spectrum.Lambda;
    assert(~isreal(lambda),'One or more eigenvalues must be underdamped for FRC computation using SSMs')
    
    % detect resonant eigenvalues in the parameter range
    switch lower(parName)
        case 'freq'
            assert(~isempty(obj.System.fext.epsilon), 'The epsilon field is empty in the dynamical system external forcing');
            [resLambda,resFreq] = find_eigs_in_freq_range(parRange,lambda);
            % obtain subintervals around each resonant eigenvalue
            [parNodes, nSubint] = subdivide_freq_range(parRange, resFreq);
        case 'amp'
            assert(~isempty(obj.System.Omega), 'The Omega field is empty in the dynamical system');
            % find eigenvalue nearest to forcing frequency
            [~,idx] = min(abs(lambda - 1i*obj.System.Omega));
            resLambda = lambda(idx);
            % setup for loop
            parNodes = [parRange(1) parRange(end)];
            nSubint = 1;
    end
    
    FRC = cell(nSubint,1);
    for i=1:nSubint
        % tune subinterval
        parSubRange = parNodes(i:i+1)';
        parSubRange = tune_parameter_range(parSubRange, obj.FRCOptions.frac, i, nSubint);
        
        % detect modes resonant with resLambda(i)
        [resModes,mFreqs] = detect_resonant_modes(resLambda(i),lambda, obj.Options.IRtol);
        
        %% FRC computation within the subinterval
        disp('*****************************************');
        disp(['Calculating FRC using SSM with master subspace: [' num2str(resModes(:).') ']']);
        
        switch obj.FRCOptions.method
            case 'level set'
                FRC{i} = obj.FRC_level_set(resModes,order,parName,parSubRange);
                plotStyle = 'circles';
            case 'continuation ep'
                % call continuation based method
                mFreqs = mFreqs(1:2:end)';
                runid  = ['freqSubint',num2str(i)];
                FRC{i}   = obj.FRC_cont_ep(runid,resModes,order,mFreqs,parName,parSubRange);
                plotStyle = 'lines';
            case 'continuation po'
                error('not implemented yet')
        end
    end
    % concatenate cell contents as struct arrays
    FRC = cat(1,FRC{:});
    
    totalComputationTime(j) = toc(startFRC);
    
    %% plot FRC in physical coordinates
    plot_FRC_full(FRC,obj.FRCOptions.outdof,order,parName,plotStyle, figs, colors(j,:))
    
end


for j = 1:numel(ORDER)    
    disp(['Total time spent on FRC computation upto O(' num2str(ORDER(j)) ') = ' datestr(datenum(0,0,0,0,0,totalComputationTime(j)),'HH:MM:SS')])
end

end

function parRange = tune_parameter_range(parRange, frac, i, nSubint)
% amplify the parameter subinterval except on the first and the last nodes.
if i>1
    parRange(1) = frac(1)*parRange(1);
end
if i<nSubint
    parRange(2) = frac(2)*parRange(2);
end
end

function [resLambda, resFreq] = find_eigs_in_freq_range(Omega,lambda)
natFreq = imag(lambda);
resFreqID = intersect(find(natFreq>Omega(1)), find(natFreq<Omega(end)));
resFreq = natFreq(resFreqID);
resLambda = lambda(resFreqID);
assert(~isempty(resFreq),'Input frequency range should include at least one natural frequency'); % we could still program this case
end

function [freqNodes, nSubint] = subdivide_freq_range(parRange,resFreq)
nSubint = numel(resFreq);
freqNodes = 0.5*(resFreq(1:end-1)+resFreq(2:end)); % center points of two adjacent resonant modes
freqNodes = [parRange(1); freqNodes; parRange(2)];
end