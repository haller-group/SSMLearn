function [NInfo,damps,freqs,str_eqn] = polarNormalForm(NInfo,optPlot,varargin) 
% Compute instantanueous damping and frequency as well as display the
% normal form.
%
% The amplitudes are denoted \rho\in\mathbb{R}^{k/2}, while angles are
% \theta\in\mathbb{R}^{k/2} for which z = \rho.*exp(1i*\theta)
% The complex valued normal form dynamics is transformed into polar
% coordinates so that the damping rates a(\rho) & frequencies w(\rho)
% are defined as (Dt is the sampling time)
%
% flow:    damps(\rho,\theta) = \dot{\rho}./rho,
%          freqs(\rho,\theta) = \dot{\theta}
% map :    damps(\rho,\theta) = 1/Dt\log(\rho_{k+1}/\rho_k),
%          damps(\rho,\theta) = (\theta_{k+1}-\theta_{k})/Dt
%
% Note that for nonresonant eigenvalues, damping and frequencies just
% depend on the amplitudes not on the frequencies
% If varargin is empty, the function assumes to deal with a flow.
% If varargin is nonempty, the function assumes to deal with a map
% and varargin is set to be the sampling time Dt.

coeff = NInfo.coefficients;
exponents = NInfo.exponents;
phi = NInfo.phi;

% Dynamics in terms of amplitude and frequencies
ndof = size(coeff,1); r = sym('r',[1 ndof],'real'); pos_coeff = abs(coeff)>0;
amp_exponents = exponents-transpose([pos_coeff; zeros(size(coeff))]);
disc_phase = sum(sum(abs(amp_exponents(:,1:ndof)-amp_exponents(:,ndof+1:end))));
if disc_phase>0
    N = @(r,th) coeff*phi([r.*exp(1i*th);r.*exp(-1i*th)]);
    if isempty(varargin)==1 % Flow case
        damps = @(r,th) real(N(r,th).*exp(-1i*th)./r);
        freqs = @(r,th) imag(N(r,th).*exp(-1i*th)./r);
    else
        Dt = varargin{:}; % Map case
        damps = @(r,th) log(abs(N(r,th))./r)/Dt;
        freqs = @(r,th) angle(N(r,th)./r)/Dt;
    end
else
    amp_exponents = 2*amp_exponents(:,1:ndof);
    r_pow = prod(r.^amp_exponents,2);
    aw = coeff*r_pow;
    if isempty(varargin)==1 % Flow case
        a_sym = real(aw); w_sym = imag(aw);
        damps = matlabFunction(a_sym,'Vars', {transpose(r)});
        freqs = matlabFunction(w_sym,'Vars', {transpose(r)});
    else % Map case
        Dt = varargin{:};
        a_sym = log(abs(aw))/Dt; w_sym = angle(aw)/Dt;
        damps = matlabFunction(a_sym,'Vars', {transpose(r)});
        freqs = matlabFunction(w_sym,'Vars', {transpose(r)});
    end
end

% Plot normal form in latex style on matlab figure (only for flows)
if isempty(varargin)==1 % Flow case
    if optPlot==1
        disp('Plotting figure with the polar normal form equations ...')
    end
    % Transform coeffs for amps and phase
    amp_exponents = exponents(:,1:ndof)+exponents(:,ndof+1:end);
    phase_coefficients = exponents(:,1:ndof)-exponents(:,ndof+1:end);
    terms_per_row = 3;
    max_terms = max(sum(abs(coeff)>0,2));
    terms_cell = cell(max_terms,2*ndof);
    n_amps = zeros(ndof,1); n_terms = zeros(ndof,1);
    for ii = 1:ndof
        idx_vec = zeros(1,ndof); idx_vec(ii) = 1;
        amp_exponents_i = amp_exponents - ...
            repmat(idx_vec,size(coeff,2),1);
        phase_coefficients_i = phase_coefficients - ...
            repmat(idx_vec,size(coeff,2),1); cc = 0;
        coeff_class  = []; % 1:amp only, 2:phase
        for jj = 1:size(coeff,2)
            if abs(coeff(ii,jj))>0
                cc = cc + 1;
                % Get the term
                if sum(abs(phase_coefficients_i(jj,:))>0)==0 % Only amp. dependence
                    coeff_class  = [coeff_class; 1];
                    monomial = monomialPNF(amp_exponents_i(jj,:));
                    if sign(real(coeff(ii,jj))) == 1
                        terms_cell{cc,ii} = ['+' num2str(abs(real(coeff(ii,jj)))) monomial];
                    else
                        terms_cell{cc,ii} = ['-' num2str(abs(real(coeff(ii,jj)))) monomial];
                    end
                    if sign(imag(coeff(ii,jj))) == 1
                        terms_cell{cc,ii+ndof} = ['+' num2str(abs(imag(coeff(ii,jj)))) monomial];
                    else
                        terms_cell{cc,ii+ndof} = ['-' num2str(abs(imag(coeff(ii,jj)))) monomial];
                    end
                else
                    coeff_class  = [coeff_class; 2];
                    monomial = monomialPNF(amp_exponents_i(jj,:),phase_coefficients_i(jj,:));
                    terms_cell{cc,ii} = ['+\mathrm{Re}((' num2str(coeff(ii,jj)) ')' monomial ')'];
                    terms_cell{cc,ii+ndof} = ['+\mathrm{Im}((' num2str(coeff(ii,jj)) ')' monomial ')'];
                end
            end
        end
        n_amps(ii) = sum(coeff_class==1); n_terms(ii) = length(coeff_class);
        [~,pos] = sort(coeff_class);
        terms_cell(:,ii) = terms_cell(pos,ii);
        terms_cell(:,ii+ndof) = terms_cell(pos,ii+ndof);
    end
    % Plot
    disc_table = mod(max_terms,terms_per_row);
    if disc_table == 0
        nrows_per_eq = max_terms/terms_per_row;
    else
        nrows_per_eq = floor(max_terms/terms_per_row)+1;
    end
    cell_str_eqn = cell(1,2*ndof); str_eqn = [];
    if ndof ==1
        str_eqn = [str_eqn '\dot{\rho}\rho^{-1} = '];
        for jj = 1:n_terms
            if jj<(n_amps+1)
                if mod(jj,terms_per_row) == 1
                    if jj == 1
                        str_eqn = [str_eqn '&'];
                    else
                        str_eqn = [str_eqn '\\ &'];
                    end
                end
                str_eqn = [str_eqn terms_cell{jj,1}];
            else
                str_eqn = [str_eqn '\\ &' terms_cell{jj,1}];
            end
        end
        cell_str_eqn{1} = str_eqn; str_eqn = [];
        str_eqn = [str_eqn '\\ \dot{\theta} = '];
        for jj = 1:n_terms
            if jj<(n_amps+1)
                if mod(jj,terms_per_row) == 1
                    if jj == 1
                        str_eqn = [str_eqn '&'];
                    else
                        str_eqn = [str_eqn '\\ &'];
                    end
                end
                str_eqn = [str_eqn terms_cell{jj,2}];
            else
                str_eqn = [str_eqn '\\ &' terms_cell{jj,2}];
            end
        end
        cell_str_eqn{2} = str_eqn; str_eqn = [];
    else
        for ii = 1:ndof
            if ii ==1
                str_eqn = [str_eqn '\dot{\rho}_{' num2str(ii) '}\rho^{-1}_{' num2str(ii) '} = '];
            else
                str_eqn = [str_eqn '\\ \dot{\rho}_{' num2str(ii) '}\rho^{-1}_{' num2str(ii) '} = '];
            end
            for jj = 1:n_terms(ii)
                if jj<(n_amps(ii)+1)
                    if mod(jj,terms_per_row) == 1
                        if jj == 1
                            str_eqn = [str_eqn '&'];
                        else
                            str_eqn = [str_eqn '\\ &'];
                        end
                    end
                    str_eqn = [str_eqn terms_cell{jj,ii}];
                else
                    str_eqn = [str_eqn '\\ &' terms_cell{jj,ii}];
                end
            end
            cell_str_eqn{ii} = str_eqn; str_eqn = [];
        end
        for ii = 1:ndof
            str_eqn = [str_eqn '\\ \dot{\theta}_{' num2str(ii) '} = '];
            for jj = 1:n_terms(ii)
                if jj<(n_amps(ii)+1)
                    if mod(jj,terms_per_row) == 1
                        if jj == 1
                            str_eqn = [str_eqn '&'];
                        else
                            str_eqn = [str_eqn '\\ &'];
                        end
                    end
                    str_eqn = [str_eqn terms_cell{jj,ii+ndof}];
                else
                    str_eqn = [str_eqn '\\ &' terms_cell{jj,ii+ndof}];
                end
            end
            cell_str_eqn{ii+ndof} = str_eqn; str_eqn = [];
        end
    end
    str_eqn = [str_eqn '\end{array}$'];
    
    switch ndof
        case 1
            str_not = '$z = \rho e^{i\theta}$';
        case 2
            str_not = '$z_1 = \rho_1 e^{i\theta_1}, \, z_2 = \rho_2 e^{i\theta_2}$';
        otherwise
            str_not = ['$z\in\mathbf{C}^' num2str(ndof) ', \,\, z_j = \rho_j e^{i\theta_j}\,\,$ for $j = 1,...,\,' num2str(ndof) '$'];
            
    end
    s_begin = '\begin{array}{rl}'; s_end = '\end{array}';
    str_eqn = [s_begin cell_str_eqn{:} s_end];
    if length(str_eqn)<1133 % Character limit is 1137
        str_eqn_plot = str_eqn; fprintf('\b Done. \n')
    else
        str_eqn_plot = [s_begin cell_str_eqn{1} s_end];
        s_len = length(str_eqn_plot);
        if  s_len<1125
            for ii = 2:2*ndof
                if (s_len + length(cell_str_eqn{ii}))<1125
                    str_eqn_plot = [s_begin cell_str_eqn{1:ii} s_end];
                    s_len = length(str_eqn_plot);
                else
                    disp(['System of equations too large for complete plotting. The figure displays the first ' num2str(ii-1) ' equations; see the variable for the full system.'])
                    str_eqn_plot = [s_begin cell_str_eqn{1:(ii-1)} '\\ & ...' s_end];
                    break;
                end
            end
        else
            str_eqn_plot = 0;
        end
    end
    if optPlot==1
        if length(str_eqn_plot)>1
            figure;
            h = plot(0,0);
            set(gcf,'color','w');
            str_above = ['Using the notation ' str_not ', the normal form reads'];
            annotation('textbox','FontSize',18,'Interpreter','latex','FaceAlpha','1','EdgeColor','w','Position',[0.01 0.1 0.99 0.9], 'String',str_above);
            annotation('textbox','FontSize',18,'Interpreter','latex','FaceAlpha','1','EdgeColor','w','Position',[0.02 0.12 0.98 0.76],'String',['$' str_eqn_plot '$']);
            delete(h);
            set(gca,'Visible','off')
        else
            disp('System of equations too large for complete plotting. See the variable for the full system.')
        end
    end
    NInfo.LaTeXPolar = str_eqn;
end
NInfo.damping = damps;
NInfo.frequency = freqs;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function monomial = monomialPNF(amp_exponents,varargin)
ndof = length(amp_exponents); monomial = [];
if isempty(varargin) == 1
    if ndof == 1
        monomial = [monomial rktostr(' \rho',amp_exponents,[])];
    else
        for ii = 1:ndof
            monomial = [monomial rktostr(' \rho',amp_exponents(ii),['_{' num2str(ii) '}'])];
        end
    end
else
    phase_coefficients = varargin{:};
    if ndof == 1
        monomial = [monomial rktostr(' \rho',amp_exponents,[])];
        monomial = [monomial ' e^{i(' tktostr(' \theta',phase_coefficients,[]) ')}'];
    else
        for ii = 1:ndof
            monomial = [monomial rktostr(' \rho',amp_exponents(ii),['_{' num2str(ii) '}'])];
        end
        monomial = [monomial ' e^{i('];
        for ii = 1:ndof
            if sign(phase_coefficients(ii)) == 1
                monomial = [monomial '+' tktostr(' \theta',abs(phase_coefficients(ii)),['_{' num2str(ii) '}'])];
            else
                monomial = [monomial '-' tktostr(' \theta',abs(phase_coefficients(ii)),['_{' num2str(ii) '}'])];
            end
        end
        monomial = [monomial ')}'];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function z_str = rktostr(z_name,exponent,suff)
switch exponent
    case 0
        z_str = '';
    case 1
        z_str = [z_name suff];
    otherwise
        z_str = [z_name '^{' num2str(exponent) '}' suff];
end
end

function z_str = tktostr(z_name,exponent,suff)
switch exponent
    case 0
        z_str = '';
    case 1
        z_str = [z_name suff];
    otherwise
        z_str = [num2str(exponent) z_name  suff];
end
end


