function [I_matrix] = I_matrix_fn(p)

% returns matrix for light integration
I_matrix = zeros((p.Xn-1)*(p.Zn-1));

%% old implementation
if(false)
    for row=1:p.Zn-1 % row nr
        
        % the top row, all columns in matrix
        index_x  = 1 +(row-1)*(p.Xn-1);
        index_z  = (1:p.Xn-1) +(row-1)*(p.Xn-1);
        I_matrix(index_z, index_x) = p.Z_vol(1,:)';
        
        for j = 2:row % loop over all the rows in the matrix corresponding to row nr: row
            index_x  = j +(row-1)*(p.Xn-1);
            index_z  = (1:p.Xn-1) +(row-1)*(p.Xn-1);
            I_matrix(index_z, index_x) = p.Z_vol(j,:)' - p.Z_vol(j-1,:)';
        end
    end
end

%% New implementation
% better approximation of the light attenuation integral
if(true)
    for row=1:p.Zn-1 % row nr
        for j = 1:row-1 % loop over all the rows in the matrix corresponding to row nr: row
            index_x  = j +(row-1)*(p.Xn-1);
            index_z  = (1:p.Xn-1) +(row-1)*(p.Xn-1);
            I_matrix(index_z, index_x) = 0.5*(p.Z(j+1,(1:end-1)) + p.Z(j+1,(2:end)))' - 0.5*(p.Z(j,(1:end-1)) + p.Z(j,(2:end)))';
        end 
        % row correspoding to depth being integrated to 
        index_x  = row +(row-1)*(p.Xn-1);
        index_z  = (1:p.Xn-1) +(row-1)*(p.Xn-1);
        I_matrix(index_z, index_x) = p.Z_vol(row,:)' - 0.5*(p.Z(row,(1:end-1)) + p.Z(row,(2:end)))';
    end
end

end

