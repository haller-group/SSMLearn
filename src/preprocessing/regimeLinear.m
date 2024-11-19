function index_time_linear = regimeLinear(data,lim)

index_linear_traj = zeros(size(data,1),1);
for jj = 1:size(data,1)
    t = data{jj,1};
    x = data{jj,2};

    index_linear_start_vect = zeros(size(x,1),1);
    for ii = 1:size(x,1)
        [amp,freq,~,time] = PFFk(t,x(ii,:),1); 
        signal_freq = freq(~isnan(freq));
        signal_amp = amp(~isnan(freq));
        signal_time = time(~isnan(freq));
        [index_linear_start_curr, ~] = fun_linear_regime(signal_freq,lim);

        index_linear_start_vect(ii) = index_linear_start_curr;
    end    
    index_linear_start_curr = round(mean(index_linear_start_vect));
    signal_time_linear = signal_time(index_linear_start_curr:end);
    [~,index_time_linear_curr] = min(abs(t-signal_time_linear(1)));

    index_linear_traj(jj) = round(mean(index_time_linear_curr));
end
index_time_linear = round(mean(index_linear_traj));

end

function [index_linear_start, lim_value] = fun_linear_regime(signal,lim)
    lim_value = lim;
    jj = 0;
    mean_vect = zeros(1,length(signal));
    for ii = 1:length(signal)-1
        mean_curr = mean(signal(1+jj:end));
        mean_vect(ii) = mean_curr;
        jj = jj + 1;
    end
    mean_vect(end) = signal(end);
    vel_mean_vect = diff(mean_vect);

    vel_abs_mean_vect = abs(vel_mean_vect);
        [~,index_linear] = find(vel_abs_mean_vect<lim_value);

    while isempty(find(vel_abs_mean_vect<lim_value))
        lim_value = exp(log(lim_value)+1);
    end
    [~,index_linear] = find(vel_abs_mean_vect<lim_value);
    index_linear_start = index_linear(1);
end
