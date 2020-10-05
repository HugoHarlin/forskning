function [dAdt] = rhs_circmodel(t,A, p)
% rhs of the circ model and standard diffusion

A = reshape(A,p.n,p.n);
dAdt = zeros(p.n);

for i=1:p.n % row nr
    for j=1:p.n % kolumn nr
        
        %%%%%%%%%%%%%%% Circulation model %%%%%%%%%%%%%%%%%5
        
        temp = 0;
        
        if((mod(i,2) == 0 && mod(j,2) == 0) || (mod(i,2)~= 0 && mod(j,2) ~= 0)) % both row and column nr odd or even, node type 1
            
            
            if(j>1) % not on the left border
                temp = temp -A(i,j); % flow to the left
            end
            
            if(j<p.n) % not on the right border
                temp = temp -A(i,j); % flow to the right
            end
            
            if(i>1) % not on the top row
                temp = temp + A(i-1,j); % inflow from the node above
            end
            
            if(i<p.n) % not on the bottom row
                temp = temp + A(i+1,j); % inflow from the node below
            end
            

        else % column nr odd and row nr even or vice versa, node type 2
            
            if(j>1) % not on the left border
                temp = temp + A(i,j-1); % inflow from the left
            end
            
            if(j<p.n) % not on the right border
                temp = temp + A(i,j+1); % inflow to the right
            end
            
            if(i>1) % not on the top row
                temp = temp - A(i,j); % flow to the node above
            end
            
            if(i<p.n) % not on the bottom row
                temp = temp - A(i,j); % blow to the node below
            end
        end
                
        %dAdt(i,j) = p.v*p.c*temp; % rate of change of A in node (i,j). c*v = diff_coeff
        dAdt(i,j) = 2*p.diff_coeff/(p.c*p.L)*temp; % rate of change of A in node (i,j). c*v = diff_coeff        
        %dAdt(i,j) = p.v/p.L*p.c*temp; % rate of change of A in node (i,j). c*v = diff_coeff        
        
        
    end
end

dAdt = dAdt(:);
end


