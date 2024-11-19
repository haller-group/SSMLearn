function [data_projected, P_min, V_trunc_slow,data_non_projected,data_non_projected_trunc] = obliqueProjection(data,index_proj,SSMDim, overEmbed, ShiftStep, varargin)
    
    time_mat = cell2mat(data(:,1));
    data_mat = cell2mat(data(:,2));

    flag_end = varargin{1,2};
    oblique_projection = varargin{1,4};

    % truncation of the final part of the data
    if flag_end == 0
        index_end = length(time_mat(1,:));
    else
        index_end = flag_end * index_proj;
    end

    time_mat_trunc = time_mat(:,index_proj:index_end);   
    data_mat_trunc = data_mat(:,index_proj:index_end);
    
    data_trunc = cell(size(data));
    m = size(time_mat_trunc,1);
    p = size(data{1,2},1);
    mp_times = ones(1,m);
    data_trunc(:,1) = mat2cell(time_mat_trunc,mp_times,size(time_mat_trunc,2));
    mp = p * ones(1,m);
    data_trunc(:,2) = mat2cell(data_mat_trunc,mp,size(time_mat_trunc,2));

    % computation of the linear oblique projection via minimization of the
    % oscillations of the backbone curve 
   
    [yData_delay, ~]=coordinatesEmbedding(data_trunc,SSMDim,'OverEmbedding',overEmbed,'ShiftSteps',ShiftStep);
    
    N = size(yData_delay{1,2},1);

    % compute the tangent space 
    dt = mean(mean(diff(time_mat_trunc,[],2)));
    X = [];
    Y = [];
    for ii = 1:1 %size(yData_delay,1)
        X = [X [yData_delay{ii,2}(:,1:end-1)]]; 
        Y = [Y [yData_delay{ii,2}(:,2:end)]];
    end
    [U,S,V] = svds(X,4);
    U = U(:,1:2*SSMDim);
    S = S(1:2*SSMDim,1:2*SSMDim);
    V = V(:,1:2*SSMDim);
    S_tilde = U.'*Y*V*pinv(S);
    [E,~,Lambda] = eigSorted(S_tilde);
    lambda_A_from_S_tilde = 1/dt*log(Lambda);  
    [~,index_sort] = sort(abs(real(lambda_A_from_S_tilde)));
    lambda_trunc = lambda_A_from_S_tilde(index_sort);
    E_sorted = E(:,index_sort);
    V_trunc = U * E_sorted;
    V_trunc_slow = [real(V_trunc(:,1)) imag(V_trunc(:,1))];

    if ~oblique_projection
        P_min = eye(N);
    else
        % minimization of the oscillations 
        var0 = compute_var(yData_delay);
        omega_0 = imag(lambda_trunc(1));
        J = @(B) minimize_oscillations_bc_variance(V_trunc_slow,B,yData_delay,var0);
        options = optimoptions('fminunc','Display','iter');
        B_0 = V_trunc_slow;
        [B_min,~] = fminunc(J,B_0,options);

        % linear oblique projection 
        P_min = V_trunc_slow*((B_min.'*V_trunc_slow)\B_min.');
    end
    
    [data_non_projected, ~]= coordinatesEmbedding(data,SSMDim,'OverEmbedding',overEmbed,'ShiftSteps',ShiftStep);
    data_projected = data_non_projected;
    for ii = 1:size(data,1)
        data_projected{ii,2} = P_min * data_non_projected{ii,2};
    end
    data_non_projected_trunc = yData_delay;


end


%% additional functions 

function [var0, var0_vect] = compute_var(data)

var0_vect = [];
lim_pff = 1e-4;
for ii = 1:size(data{1,2},1)
    kmean = 1;
    [amp,freq,~,~] = PFFk(data{1,1},data{1,2}(ii,:),kmean);
    amp = amp(~isnan(freq));
    signal = freq(~isnan(freq));
    amp = amp(amp>lim_pff);
    signal = signal(amp>lim_pff);
    average = mean(signal);
    var0_curr = var(signal - average);
    var0_vect = [var0_vect var0_curr];
end
var0 = max(var0_vect);

end

function J = minimize_oscillations_bc_variance(V_slow,B,data,var0)
    J_vect = zeros(size(data,1),1);
    for ii = size(data,1)
        sol = data{ii,2};
        P = V_slow*((B.'*V_slow)\B.');
        z = P*sol;
        data_proj = cell(1,2);
        data_proj{1,1} = data{ii,1};
        data_proj{1,2} = z;
        [~, var_vect] = compute_var(data_proj);
        J_vect(ii) = sum(var_vect(~isnan(var_vect)))/var0;
    end
    J = mean(J_vect);
end

