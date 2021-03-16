function table_nf = disp_normalform(coeffs,exponents,varargin)
% Display the obtained normal form for oscillators. The monomial refers to
% powers of complex numbers z and their conjugated z`

cell_table = cell(size(coeffs,1)+1,size(coeffs,2)+1);
cell_table{1,1} = 'Monomial:';
if size(exponents,2)==2
    if isempty(varargin)==1
        cell_table{2,1} = '\dot{z} =';
    else
        cell_table{2,1} = 'z <--|';
    end
    for ii = 1:size(exponents,1)
        cell_table{2,ii+1} = coeffs(ii);
        cell_table{1,ii+1} = [zktostr('z',exponents(ii,1),[]) ...
            zktostr('z`',exponents(ii,2),[])];
    end
else
    if isempty(varargin)==1
        for jj = 1:size(coeffs,1)
            cell_table{jj+1,1} = ['\dot{z}_' num2str(jj) ' ='];
        end
    else
        for jj = 1:size(coeffs,1)
            cell_table{jj+1,1} = ['z_' num2str(jj) ' <--|'];
        end
    end
    for ii = 1:size(exponents,1)
        name_i = [];
        for jj = 1:size(coeffs,1)
            cell_table{jj+1,1+ii} = coeffs(jj,ii);
            name_i = [name_i zktostr('z',exponents(ii,jj),['_' num2str(jj)])...
                zktostr('z`',exponents(ii,jj+size(coeffs,1)), ...
                ['_' num2str(jj)])];
        end
        cell_table{1,ii+1} = name_i ;
    end
end

table_nf = cell2table(cell_table(2:end,2:end),'VariableNames',...
                       cell_table(1,2:end),'RowNames',cell_table(2:end,1));
end

function z_str = zktostr(z_name,exponent,suff)
switch exponent
    case 0
        z_str = '';
    case 1
        z_str = [z_name suff];
    otherwise
        z_str = [z_name num2str(exponent) suff];
end
end