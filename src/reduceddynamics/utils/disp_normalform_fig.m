function [str_eqn,flag_long] = disp_normalform_fig(coeffs,exponents,termsprow,varargin)
% Display the obtained normal form for oscillators. The monomial refers to
% powers of complex numbers z and their conjugated z`. If varargin is empty
% then the function assumes a flow, otherwise a map.

% Plot normal form in latex style on matlab figure
ndof = size(coeffs,1);
terms_per_row = termsprow;
max_terms = max(sum(abs(coeffs)>0,2));
terms_cell = cell(max_terms,ndof);
for ii = 1:ndof
    cc = 0;
    for jj = 1:size(coeffs,2)
        if abs(coeffs(ii,jj))>0
            cc = cc + 1;
            monomial = monomial_nf(ndof,exponents(jj,:));
            terms_cell{cc,ii} = ['(' num2str(coeffs(ii,jj)) ')' monomial ];
            terms_cell{cc,ii+ndof} = ['(' num2str(coeffs(ii,jj)) ')' monomial ];
        end
    end
end
% Plot
disc_table = mod(max_terms,terms_per_row);
if disc_table == 0
    nrows_per_eq = max_terms/terms_per_row;
else
    nrows_per_eq = floor(max_terms/terms_per_row)+1;
end
cell_str_eqn = cell(1,ndof); str_eqn = [];
if ndof ==1
    if isempty(varargin) == 1
        str_eqn = [str_eqn '\dot{z} = '];
    else
        str_eqn = [str_eqn 'z \leftarrow '];
    end
    for jj = 1:max_terms
        if terms_per_row>1
            if jj ==1
                str_eqn = [str_eqn '&' terms_cell{jj,1}];
            else
                if mod(jj,terms_per_row) == 1
                    str_eqn = [str_eqn '+ \\ &' terms_cell{jj,1}];
                else
                    str_eqn = [str_eqn ' + ' terms_cell{jj,1}];
                end
            end
        else
            if jj == 1
                str_eqn = [str_eqn '&' terms_cell{jj,1}];
            else
                str_eqn = [str_eqn '+ \\ &' terms_cell{jj,1}];
            end
        end
    end
    cell_str_eqn{1} = str_eqn; str_eqn = [];
else
    for ii = 1:ndof
        if isempty(varargin) == 1
            if ii == 1
                str_eqn = [str_eqn '\dot{z}_{' num2str(ii) '} = '];
            else
                str_eqn = [str_eqn '\\ \dot{z}_{' num2str(ii) '} = '];
            end
        else
            if ii == 1
                str_eqn = [str_eqn 'z_{' num2str(ii) '} \leftarrow '];
            else
                str_eqn = [str_eqn '\\ z_{' num2str(ii) '} \leftarrow '];
            end
        end
        for jj = 1:max_terms
            if terms_per_row>1
                if jj ==1
                    str_eqn = [str_eqn '&' terms_cell{jj,ii}];
                else
                    if mod(jj,terms_per_row) == 1
                        str_eqn = [str_eqn '+ \\ &' terms_cell{jj,ii}];
                    else
                        str_eqn = [str_eqn ' + ' terms_cell{jj,ii}];
                    end
                end
                
            else
                if jj == 1
                    str_eqn = [str_eqn '&' terms_cell{jj,ii}];
                else
                    str_eqn = [str_eqn '+ \\ &' terms_cell{jj,ii}];
                end
            end
        end
        cell_str_eqn{ii} = str_eqn; str_eqn = [];
    end
end
s_begin = '\begin{array}{rl}'; s_end = '\end{array}';
str_eqn = [s_begin cell_str_eqn{:} s_end]; flag_long = 1;
if length(str_eqn)<1133 % Character limit is 1137
    str_eqn_plot = str_eqn; flag_long = 0;
else
    str_eqn_plot = [s_begin cell_str_eqn{1} s_end];
    s_len = length(str_eqn_plot);
    if  s_len<1125
        for ii = 2:2*ndof
            if (s_len + length(cell_str_eqn{ii}))<1125
                str_eqn_plot = [s_begin cell_str_eqn{1:ii} s_end];
                s_len = length(str_eqn_plot);
            else
                str_eqn_plot = [s_begin cell_str_eqn{1:(ii-1)} '\\ & ...' s_end];
                break;
            end
        end
    else
        str_eqn_plot = 0;
    end
end
if length(str_eqn_plot)>1
    figure;
    set(gcf,'color','w');
    str_above = ['Using the notation $\bar{\,}$ for the complex conjugated, the identified normal form is'];
    annotation('textbox','FontSize',18,'Interpreter','latex','FaceAlpha','1','EdgeColor','w','Position',[0.01 0.1 0.99 0.9], 'String',str_above);
    annotation('textbox','FontSize',18,'Interpreter','latex','FaceAlpha','1','EdgeColor','w','Position',[0.02 0.12 0.98 0.76],'String',['$' str_eqn_plot '$']);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function monomial = monomial_nf(ndof,exponents)
monomial = [];
if ndof == 1
    monomial = [monomial zktostr(' z',exponents(1),[]) ...
        zktostr(' \bar{z}',exponents(2),[])];
else
    for ii = 1:ndof
        s_i = ['_{' num2str(ii) '}'];
        monomial = [monomial zktostr(' z',exponents(ii),s_i) ...
            zktostr(' \bar{z}',exponents(ii+ndof),s_i)];
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_str = zktostr(z_name,exponent,suff)
switch exponent
    case 0
        z_str = '';
    case 1
        z_str = [z_name suff];
    otherwise
        z_str = [z_name '^{' num2str(exponent) '}' suff];
end
end
