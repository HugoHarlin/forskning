function [dAdt] = rhs_diff(t,A, p)
% rhs of the circ model and standard diffusion

dA_diffdt = zeros(p.n);
dAdt = zeros(p.n);

for i=1:p.n % row nr
    for j=1:p.n % kolumn nr
        
        %%%%%%%%%%%%%%% Circulation model %%%%%%%%%%%%%%%%%5
        
        temp = 0;
        
        if((mod(i,2) == 0 && mod(j,2) == 0) || (mod(i,2)~= 0 && mod(j,2) ~= 0)) % both row and column nr odd or even, node type 1
            
            
            if(j>1) % not on the left border
                temp = temp -A(i,j); % flow to the left
            end
            
            if(j<n) % not on the right border
                temp = temp -A(i,j); % flow to the right
            end
            
            if(i>1) % not on the top row
                temp = temp + A(i-1,j); % inflow from the node above
            end
            
            if(i<n) % not on the bottom row
                temp = temp + A(i+1,j); % inflow from the node below
            end
            
            
            
        else % column nr odd and row nr even or vice versa, node type 2
            
            if(j>1) % not on the left border
                temp = temp + A(i,j-1); % inflow from the left
            end
            
            if(j<n) % not on the right border
                temp = temp + A(i,j+1); % inflow to the right
            end
            
            if(i>1) % not on the top row
                temp = temp - A(i,j); % flow to the node above
            end
            
            if(i<n) % not on the bottom row
                temp = temp - A(i,j); % blow to the node below
            end
        end
        
        dAdt(i,j) = v/L*temp; % updating the concentration in node (i,j)
        
        
        %%%%%%%%%%%%%%%%%% Standard diffusion %%%%%%%%%%%%%%%%%%
        
        if(i==1)
            % top row
            da2dz2 = A_diff(i+1,j)-A_diff(i,j);
        elseif(i==n)
            % bottom row
            da2dz2 = A_diff(i,j)-A_diff(i-1,j);
        else
            % neither top or bottom row
            da2dz2 = A_diff(i+1,j)+ A_diff(i-1,j) -2*A_diff(i,j);
        end
        
        if(j==1)
            % left column
            da2dx2 = A_diff(i,j+1)-A_diff(i,j);
            
        elseif(j==n)
            % right column
            da2dx2 = A_diff(i,j)-A_diff(i,j-1);
        else
            % neither left nor right column
            da2dx2 = A_diff(i,j+1) + A_diff(i,j-1) -2*A_diff(i,j);
        end
        
        %da2dz2 = 0;
        % da2dx2 = 0;
        dA_diffdt(i,j) = v/L*(da2dx2 + da2dz2); % updating the concentration in node (i,j)
        
    end
end

end


