function [I_matrix] = I_matrix_fn(p)
% returns matrix for light integration


I_matrix = zeros((p.Xn-1)*(p.Zn-1));

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

