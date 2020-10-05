function [dA_diffdt] = rhs_standard_diff(t,A_diff,p)
% rhs of the circ model and standard diffusion

A_diff = reshape(A_diff,p.n,p.n);
dA_diffdt = zeros(p.n);

for i=1:p.n % row nr
    for j=1:p.n % kolumn nr
        
        %%%%%%%%%%%%%%%%%% Standard diffusion %%%%%%%%%%%%%%%%%%
        
        if(i==1)
            % top row
            da2dz2 = A_diff(i+1,j)-A_diff(i,j);
        elseif(i==p.n)
            % bottom row
            da2dz2 = -A_diff(i,j) + A_diff(i-1,j);
        else
            % neither top or bottom row
            da2dz2 = A_diff(i+1,j)+ A_diff(i-1,j) -2*A_diff(i,j);
        end
        
        if(j==1)
            % left column
            da2dx2 = A_diff(i,j+1)-A_diff(i,j);
            
        elseif(j==p.n)
            % right column
            da2dx2 = -A_diff(i,j)+A_diff(i,j-1);
        else
            % neither left nor right column
            da2dx2 = A_diff(i,j+1) + A_diff(i,j-1) -2*A_diff(i,j);
        end
        
        %da2dz2 = 0;
        % da2dx2 = 0;
        dA_diffdt(i,j) = p.diff_coeff*(1/p.L^2)*(da2dx2 + da2dz2); % updating the concentration in node (i,j)
        
    end
end
dA_diffdt = dA_diffdt(:);

end


