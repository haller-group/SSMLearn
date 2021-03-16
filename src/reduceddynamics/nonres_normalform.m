function [a,w] = nonres_normalform(coeff,exponents,varargin)
% Compute damping and frequency as well as display the normal form for
% nonresonant eigenvalues in terms of amplitude

% Dynamics in terms of amplitude and frequencies
ndof = size(coeff,1); r = sym('r',[1 ndof]); pos_coeff = abs(coeff)>0;
amp_exponents = exponents-transpose([pos_coeff; zeros(size(coeff))]);
if sum(sum(abs(amp_exponents(:,1:ndof)-amp_exponents(:,ndof+1:end))))>0
    disp('Non separable normal form.')
    a = []; w = [];
else
    amp_exponents = 2*amp_exponents(:,1:ndof);
    r_pow = prod(r.^amp_exponents,2);
    aw = coeff*r_pow;
    if isempty(varargin)==1 % Flow case
        a_sym = real(aw); w_sym = imag(aw);
        a = matlabFunction(a_sym,'Vars', {transpose(r)});
        w = matlabFunction(w_sym,'Vars', {transpose(r)});
        % Resorting coefficients for sorting
        err_flag = 0 ;
        if ndof == 1
            ua_coeff = coeff; ua_exponents = amp_exponents;
        else
            idx_i = find(abs(coeff(1,:)));
            ua_coeff = coeff(1,idx_i);
            [ua_exponents,pos] = sortrows100(amp_exponents(idx_i,:));
            ua_coeff = ua_coeff(pos);
            for ii = 2:ndof
                idx_i = find(abs(coeff(ii,:)));
                ua_coeff_i = coeff(ii,idx_i);
                [ua_exponents_i,pos] = sortrows100(amp_exponents(idx_i,:));
                if size(ua_exponents,1)~=size(ua_exponents_i)
                    err_flag = 1;
                else
                    if sum(sum(abs(ua_exponents-ua_exponents_i)))>0
                        err_flag = 1;
                    else
                        ua_coeff = [ua_coeff; ua_coeff_i(pos)];
                    end
                end
            end
        end
        if err_flag == 0 % Display the normal form
            disp(['The data-driven polar normal form dynamics reads:'])
            fprintf('\n')
            table_nf = disp_nonres_normalform_flow(ua_coeff,ua_exponents);
            disp(table_nf)
            if ndof == 1
                disp(['Notation: z = r*exp(1i*\theta) where r is the '...
                    'amplitude and \theta the phase; ' ...
                    'rk is the k-th power of r.'])
            else
                
                disp(['Notation: z_j = r_j*exp(1i*\theta_j) where r_j is '...
                    'the amplitude and \theta_j the phase; ' ...
                    'rk_j is the k-th power of r_j.'])
            end
        end
        
    else % Map case
        Dt = varargin{:};
        a_sym = log(abs(aw))/Dt; w_sym = phase(aw)/Dt;
        a = matlabFunction(a_sym,'Vars', {transpose(r)});
        w = matlabFunction(w_sym,'Vars', {transpose(r)});
        
        % Resorting coefficients for sorting
        err_flag = 0 ;
        if ndof == 1
            ua_coeff = coeff; ua_exponents = amp_exponents;
        else
            idx_i = find(abs(coeff(1,:)));
            ua_coeff = coeff(1,idx_i);
            [ua_exponents,pos] = sortrows100(amp_exponents(idx_i,:));
            ua_coeff = ua_coeff(pos);
            for ii = 2:ndof
                idx_i = find(abs(coeff(ii,:)));
                ua_coeff_i = coeff(ii,idx_i);
                [ua_exponents_i,pos] = sortrows100(amp_exponents(idx_i,:));
                if size(ua_exponents,1)~=size(ua_exponents_i)
                    err_flag = 1;
                else
                    if sum(sum(abs(ua_exponents-ua_exponents_i)))>0
                        err_flag = 1;
                    else
                        ua_coeff = [ua_coeff; ua_coeff_i(pos)];
                    end
                end
            end
        end
        if err_flag == 0 % Display the normal form
            disp(['The data-driven polar normal form dynamics reads:'])
            fprintf('\n')
            table_nf = disp_nonres_normalform_map(ua_coeff,ua_exponents);
            disp(table_nf)
            if ndof == 1
                disp(['Notation: z = r*exp(1i*\theta) where r is the '...
                    'amplitude and \theta the phase; ' ...
                    'rk is the k-th power of r.'])
                disp(['The polar form map '...
                    ' consists of r <--| |R(r)|r and ' ...
                    '\theta <--| \theta  + \arg(R(r))'])
            else
                
                disp(['Notation: z_j = r_j*exp(1i*\theta_j) where r_j is '...
                    'the amplitude and \theta_j the phase; ' ...
                    'rk_j is the k-th power of r_j.'])
                disp(['The polar form '...
                    'map consists of r_j <--| |R_j(r)|r_j and ' ...
                    '\theta_j <--| \theta_j  + \arg(R_j(r))'])
            end
        end
        
    end
end
end

function table_nf = disp_nonres_normalform_flow(coeffs,exponents)
% Disply the obtained normal form flows for non resonant oscillators. The
% monomial refers to powers of the amplitudes r

cell_table = cell(2*size(coeffs,1)+1,size(coeffs,2)+1);
cell_table{1,1} = 'Monomial:'; ndof = size(exponents,2);
if ndof==1
    cell_table{2,1} = '\dot{r} =';
    cell_table{3,1} = '\dot{\theta} =';
    for ii = 1:size(exponents,1)
        cell_table{2,ii+1} = real(coeffs(ii));
        cell_table{3,ii+1} = imag(coeffs(ii));
        if exponents(ii) == 0
            cell_table{1,ii+1} = '1';
        else
            cell_table{1,ii+1} = rktostr('r',exponents(ii),[]);
        end
    end
else
    for jj = 1:size(coeffs,1)
        cell_table{jj+1,1} = ['\dot{r}_' num2str(jj) ' ='];
        cell_table{jj+1+ndof,1} = ['\dot{\theta}_' num2str(jj) ' ='];
    end
    for ii = 1:size(exponents,1)
        if sum(exponents(ii,:))==0
            name_i = '1';
            for jj = 1:ndof
                cell_table{jj+1,ii+1} = real(coeffs(jj,ii));
                cell_table{jj+1+ndof,ii+1} = imag(coeffs(jj,ii));
            end
        else
            name_i = [];
            for jj = 1:ndof
                cell_table{jj+1,ii+1} = real(coeffs(jj,ii));
                cell_table{jj+1+ndof,ii+1} = imag(coeffs(jj,ii));
                name_i = [name_i ...
                    rktostr('r',exponents(ii,jj),['_' num2str(jj)])];
            end
        end
        cell_table{1,ii+1} = name_i ;
    end
end

table_nf = cell2table(cell_table(2:end,2:end),'VariableNames',...
    cell_table(1,2:end),'RowNames',cell_table(2:end,1));
end

function table_nf = disp_nonres_normalform_map(coeffs,exponents)
% Disply the obtained normal form map for non resonant oscillators. The
% monomial refers to powers of the amplitudes r

cell_table = cell(size(coeffs,1)+1,size(coeffs,2)+1);
cell_table{1,1} = 'Monomial:'; ndof = size(exponents,2);
if ndof==1
    cell_table{2,1} = 'R(r) =';
    for ii = 1:size(exponents,1)
        cell_table{2,ii+1} = coeffs(ii);
        if exponents(ii) == 0
            cell_table{1,ii+1} = '1';
        else
            cell_table{1,ii+1} = rktostr('r',exponents(ii),[]);
        end
    end
else
    for jj = 1:size(coeffs,1)
        cell_table{jj+1,1} = ['R_' num2str(jj) '(r) ='];
    end
    for ii = 1:size(exponents,1)
        if sum(exponents(ii,:))==0
            name_i = '1';
            for jj = 1:ndof
                cell_table{jj+1,ii+1} = coeffs(jj,ii);
            end
        else
            name_i = [];
            for jj = 1:ndof
                cell_table{jj+1,ii+1} = coeffs(jj,ii);
                name_i = [name_i ...
                    rktostr('r',exponents(ii,jj),['_' num2str(jj)])];
            end
        end
        cell_table{1,ii+1} = name_i ;
    end
end

table_nf = cell2table(cell_table(2:end,2:end),'VariableNames',...
    cell_table(1,2:end),'RowNames',cell_table(2:end,1));
end

function z_str = rktostr(z_name,exponent,suff)
switch exponent
    case 0
        z_str = '';
    case 1
        z_str = [z_name suff];
    otherwise
        z_str = [z_name num2str(exponent) suff];
end
end

function [A,pos] = sortrows100(A)
[~,pos] = sort(A*transpose(100.^[0:size(A,2)-1]),'ascend'); A = A(pos,:);
[~,pos2] = sort(sum(A,2),'ascend'); pos = pos(pos2); A = A(pos2,:);
end