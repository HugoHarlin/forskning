function [S] = Stiffness_matrix(p)

%% Stiffness matrix of the linearized system
% assembles the mass matrix of the linearized system with integrated
% boundary conditions.


S = zeros(2*(p.Xn-1)*(p.Zn-1) + 2*(p.Xn-1));

%% mass matrix for the algae
% assembling the Algae portion of the stiffness vector, looping row by row
for i = 1:p.Zn-1 % eta
    for j = 1:p.Xn-1 % xi
        
        x = j+(i-1)*(p.Xn-1);
        y = j+(i-1)*(p.Xn-1);
        
        %%%%%%%%% remove when done, only kept for reference %%%%%%%%%%
        % Diffusion and sinking of algae with  BC:s
        %         dAdt(i,j) = integral_A_1(i,j,p,A) + integral_A_2(i,j,p,A) + ...
        %             integral_A_3(i,j,p,A) + integral_A_4(i,j,p,A);
        %
        %         %integral 1
        %         output = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        %             (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi - ...
        %             p.dZ_dXi_preCalc(i,j,1) * dA_dEta_1(i,j,p,A));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Integral 1
        if(true)
            if(j ~= p.Xn-1) % we are not on the right border. If we are, the BC there implies that the stiffness matrix coefficients are zero (no flux).
                
                %  implementation of 2*pi*p.W/(p.Xn-1).*j .* p.dx * 1/p.dX_dXi_preCalc(i,j,1) * (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi
                % into the mass matrix.
                % A simple forward difference is used to approximate dA_dXi,
                % dA_dXi =  A(i,j+1) - A(i,j).
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/dX_dXi(p,i-0.5,j)*dZ_dEta(p,i-0.5,j);
                
                S(x,y) = S(x,y) - scale_factor;
                S(x,y+1) = S(x,y+1) + scale_factor;
                
                scale_factor = -2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/dX_dXi(p,i-0.5,j)*dZ_dXi(p,i-0.5,j)*0.25;
                
                if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
                    % dA_dEta = 0.25*(A(i,j) + A(i,j+1) - A(i-1,j) - A(i-1,j+1));
                    
                    S(x,y) =  S(x,y) + scale_factor;
                    S(x,y+1) =  S(x,y+1)+ scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)+1) = S(x,y+1-(p.Xn-1)) - scale_factor;
                    
                elseif(i ==  1) %  top border
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i,j) - A(i,j+1));
                    
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y+(p.Xn-1)+1) = S(x,y+(p.Xn-1)+1) + scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y+1) = S(x,y+1) - scale_factor;
                    
                else % in between top and bottom border. we can use a standard quadrature
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i-1,j) - A(i-1,j+1));
                    
                    S(x,y+(p.Xn-1))   = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y+(p.Xn-1)+1) = S(x,y+(p.Xn-1)+1) + scale_factor;
                    S(x,y-(p.Xn-1))   = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)+1) = S(x,y-(p.Xn-1)+1) - scale_factor;
                    
                end
                
            end
        end
        
        %% Integral 2
        if(true)
            if(i ~=1) % There is no flow out through the top of the lake.
                %dA_dEta = A(i,j) - A(i-1,j);
                
                % Diffusion
                %  output = 2*pi*p.W/(p.Xn-1).*(j-0.5).*(p.dx * p.dZ_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)* ...
                %     dA_dXi_2(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dZ_dEta_preCalc(i,j,2))* ...
                %    dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,2)^2 + p.dz*p.dX_dXi_preCalc(i,j,2)^2));
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.dx *dZ_dXi(p,i-1,j-0.5)/dX_dXi(p,i-1,j-0.5)*0.25;
                
                if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
                    %dA_dXi = 0.25*(A(i,j) + A(i-1,j) - A(i,j-1) - A(i-1,j-1));
                    S(x,y) =  S(x,y) + scale_factor;
                    S(x,y-(p.Xn-1)) =  S(x,y-(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-(p.Xn-1)-1) - scale_factor;
                    
                elseif(j ==  1) %  left
                    %dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j) - A(i-1,j));
                    
                    S(x,y+1) =  S(x,y+1) + scale_factor;
                    S(x,y-(p.Xn-1)+1) =  S(x,y-(p.Xn-1)+1)+ scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    
                else % in between left and right border. we can use a standard quadrature
                    %dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j-1) - A(i-1,j-1));
                    
                    S(x,y+1) =  S(x,y+1) + scale_factor;
                    S(x,y-(p.Xn-1)+1) =  S(x,y-(p.Xn-1)+1)+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-1-(p.Xn-1)) - scale_factor;
                end
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*(j-0.5)*(-1)/(dX_dXi(p,i-1,j-0.5)*dZ_dEta(p,i-1,j-0.5))*(p.dx*dZ_dXi(p,i-1,j-0.5)^2 + p.dz*dX_dXi(p,i-1,j-0.5)^2);
                S(x,y) = S(x,y) + scale_factor;
                S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                
                % Sinking
                %output = output + 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,2).*0.5.*(A(i,j)+A(i-1,j));
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* dX_dXi(p,i-1,j-0.5).*0.5;
                S(x,y) = S(x,y) + scale_factor;
                S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) + scale_factor;
                
            end
        end
        
        %% Integral 3
        if(true)
            if(j ~= 1) % we are not on the left border.
                %   dA_dXi = A(i,j)-A(i,j-1);
                % output = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/p.dX_dXi_preCalc(i,j,3) * ...
                %    (p.dZ_dEta_preCalc(i,j,3) .* dA_dXi - ...
                %    p.dZ_dXi_preCalc(i,j,3) * dA_dEta_3(i,j,p,A));
                
                scale_factor = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/dX_dXi(p,i-0.5,j-1)*dZ_dEta(p,i-0.5,j-1);
                S(x,y) = S(x,y) + scale_factor;
                S(x,y-1) = S(x,y-1) - scale_factor;
                
                scale_factor =  2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/dX_dXi(p,i-0.5,j-1)*dZ_dXi(p,i-0.5,j-1)*0.25;
                
                if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
                    %dA_dEta = 0.25*(A(i,j) + A(i,j-1) - A(i-1,j) - A(i-1,j-1));
                    S(x,y) = S(x,y) + scale_factor;
                    S(x,y-1) = S(x,y-1) + scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-(p.Xn-1)-1) - scale_factor;
                    
                elseif(i ==  1) %  top border
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i,j) - A(i,j-1));
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y-1+(p.Xn-1)) = S(x,y-1+(p.Xn-1)) + scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    
                else % in between top and bottom border. we can use a standard quadrature
                    % dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i-1,j) - A(i-1,j-1));
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y-1+(p.Xn-1)) = S(x,y-1+(p.Xn-1)) + scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-1-(p.Xn-1)) = S(x,y-1-(p.Xn-1)) - scale_factor;
                end
                
            end
        end
        
        %% Integral 4
        if(true)
            if(i ~=p.Zn-1) % No diffusive flux through the bottom border.
                %dA_dEta = A(i+1,j) - A(i,j);
                %Diffusion
                % output = -2*pi*p.W/(p.Xn-1).*(j-0.5) .* (p.dx * p.dZ_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
                %     dA_dXi_4(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,4)*p.dZ_dEta_preCalc(i,j,4))* ...
                %     dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,4)^2 + p.dz*p.dX_dXi_preCalc(i,j,4)^2));
                
                
                scale_factor =  -2*pi*p.W/(p.Xn-1).*(j-0.5).*p.dx*dZ_dXi(p,i,j-0.5)/dX_dXi(p,i,j-0.5)*0.25;
                
                if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
                    %  dA_dXi = 0.25*(A(i,j) + A(i+1,j) - A(i,j-1) - A(i+1,j-1));
                    
                    S(x,y) = S(x,y) + scale_factor;
                    S(x,y+(p.Xn-1)) =  S(x,y+(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y+(p.Xn-1)-1) = S(x,y+(p.Xn-1)-1) - scale_factor;
                    
                elseif(j ==  1) %  left
                    % dA_dXi = 0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j) - A(i+1,j));
                    S(x,y+1) = S(x,y+1) + scale_factor;
                    S(x,y+1+(p.Xn-1)) =  S(x,y+1+(p.Xn-1))+ scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) - scale_factor;
                    
                else % in between left and right border. we can use a standard quadrature
                    %dA_dXi = 0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j-1) - A(i+1,j-1));
                    S(x,y+1) = S(x,y+1) + scale_factor;
                    S(x,y+1+(p.Xn-1)) =  S(x,y+1+(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y+(p.Xn-1)-1) = S(x,y+(p.Xn-1)-1) - scale_factor;
                end
                
                scale_factor =  2*pi*p.W/(p.Xn-1).*(j-0.5).*1/( dX_dXi(p,i,j-0.5)*dZ_dEta(p,i,j-0.5))*(p.dx*dZ_dXi(p,i,j-0.5)^2 + p.dz*dX_dXi(p,i,j-0.5)^2);
                S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                S(x,y) = S(x,y) - scale_factor;
                
                % Sinking
                % output = output - 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,4).*0.5.*(A(i+1,j)+A(i,j));
                scale_factor = - 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.*dX_dXi(p,i,j-0.5).*0.5;
                S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                S(x,y) = S(x,y) + scale_factor;
                
            else
                % algae that sink into the bottom become sedimented algae
                % output = - p.death_rate*2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.* p.dX_dXi_preCalc(i,j,4)*A(i,j);
                S(x,y) = S(x,y) - p.death_rate*2*pi*p.W/(p.Xn-1).*(j-0.5).*p.v.*dX_dXi(p,i,j-0.5);
                
            end
            
        end
        
    end
end

%% mass matrix for the dissolved nutrients
% assembling the Algae portion of the stiffness vector, looping row by row
for i = 1:p.Zn-1 % eta
    for j = 1:p.Xn-1 % xi
        
        x = j+(i-1)*(p.Xn-1) + (p.Xn-1)*(p.Zn-1);
        y = j+(i-1)*(p.Xn-1) + (p.Xn-1)*(p.Zn-1);
        
        %%%%%%%%% remove when done, only kept for reference %%%%%%%%%%
        % Diffusion and sinking of algae with  BC:s
        %         dAdt(i,j) = integral_Rd_1(i,j,p,A) + integral_Rd_2(i,j,p,A) + ...
        %             integral_Rd_3(i,j,p,A) + integral_Rd_4(i,j,p,A);
        %
        %         %integral 1
        %         output = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        %             (p.dZ_dEta_preCalc(i,j,1) .* dRd_dXi - ...
        %             p.dZ_dXi_preCalc(i,j,1) * dRd_dEta_1(i,j,p,A));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Integral 1
        if(true)
            if(j ~= p.Xn-1) % we are not on the right border. If we are, the BC there implies that the stiffness matrix coefficients are zero (no flux).
                
                %  implementation of 2*pi*p.W/(p.Xn-1).*j .* p.dx * 1/p.dX_dXi_preCalc(i,j,1) * (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi
                % into the mass matrix.
                % A simple forward difference is used to approximate dA_dXi,
                % dRd_dXi =  Rd(i,j+1) - Rd(i,j).
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/dX_dXi(p,i-0.5,j)*dZ_dEta(p,i-0.5,j);
                
                S(x,y) = S(x,y) - scale_factor;
                S(x,y+1) = S(x,y+1) + scale_factor;
                
                scale_factor = -2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/dX_dXi(p,i-0.5,j)*dZ_dXi(p,i-0.5,j)*0.25;
                
                if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
                    % dRd_dEta = 0.25*(Rd(i,j) + Rd(i,j+1) - Rd(i-1,j) - Rd(i-1,j+1));
                    
                    S(x,y) =  S(x,y) + scale_factor;
                    S(x,y+1) =  S(x,y+1)+ scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)+1) = S(x,y+1-(p.Xn-1)) - scale_factor;
                    
                elseif(i ==  1) %  top border
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i,j) - A(i,j+1));
                    
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y+(p.Xn-1)+1) = S(x,y+(p.Xn-1)+1) + scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y+1) = S(x,y+1) - scale_factor;
                    
                else % in between top and bottom border. we can use a standard quadrature
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i-1,j) - A(i-1,j+1));
                    
                    S(x,y+(p.Xn-1))   = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y+(p.Xn-1)+1) = S(x,y+(p.Xn-1)+1) + scale_factor;
                    S(x,y-(p.Xn-1))   = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)+1) = S(x,y-(p.Xn-1)+1) - scale_factor;
                    
                end
                
            end
        end
        
        %% Integral 2
        if(true)
            if(i ~=1) % There is no flow out through the top of the lake.
                %dA_dEta = A(i,j) - A(i-1,j);
                
                % Diffusion
                %  output = 2*pi*p.W/(p.Xn-1).*(j-0.5).*(p.dx * p.dZ_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)* ...
                %     dA_dXi_2(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dZ_dEta_preCalc(i,j,2))* ...
                %    dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,2)^2 + p.dz*p.dX_dXi_preCalc(i,j,2)^2));
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*(j-0.5).*p.dx *dZ_dXi(p,i-1,j-0.5)/dX_dXi(p,i-1,j-0.5)*0.25;
                
                if(j ==  p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
                    %dA_dXi = 0.25*(A(i,j) + A(i-1,j) - A(i,j-1) - A(i-1,j-1));
                    S(x,y) =  S(x,y) + scale_factor;
                    S(x,y-(p.Xn-1)) =  S(x,y-(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-(p.Xn-1)-1) - scale_factor;
                    
                elseif(j ==  1) %  left
                    %dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j) - A(i-1,j));
                    
                    S(x,y+1) =  S(x,y+1) + scale_factor;
                    S(x,y-(p.Xn-1)+1) =  S(x,y-(p.Xn-1)+1)+ scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    
                else % in between left and right border. we can use a standard quadrature
                    %dA_dXi = 0.25*(A(i,j+1) + A(i-1,j+1) - A(i,j-1) - A(i-1,j-1));
                    
                    S(x,y+1) =  S(x,y+1) + scale_factor;
                    S(x,y-(p.Xn-1)+1) =  S(x,y-(p.Xn-1)+1)+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-1-(p.Xn-1)) - scale_factor;
                end
                
                scale_factor = 2*pi*p.W/(p.Xn-1).*(j-0.5)*(-1)/(dX_dXi(p,i-1,j-0.5)*dZ_dEta(p,i-1,j-0.5))*(p.dx*dZ_dXi(p,i-1,j-0.5)^2 + p.dz*dX_dXi(p,i-1,j-0.5)^2);
                S(x,y) = S(x,y) + scale_factor;
                S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                
            end
        end
        
        %% Integral 3
        if(true)
            if(j ~= 1) % we are not on the left border.
                %   dA_dXi = A(i,j)-A(i,j-1);
                % output = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/p.dX_dXi_preCalc(i,j,3) * ...
                %    (p.dZ_dEta_preCalc(i,j,3) .* dA_dXi - ...
                %    p.dZ_dXi_preCalc(i,j,3) * dA_dEta_3(i,j,p,A));
                
                scale_factor = -2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/dX_dXi(p,i-0.5,j-1)*dZ_dEta(p,i-0.5,j-1);
                S(x,y) = S(x,y) + scale_factor;
                S(x,y-1) = S(x,y-1) - scale_factor;
                
                scale_factor =  2*pi*p.W/(p.Xn-1).*(j-1) .*  p.dx * 1/dX_dXi(p,i-0.5,j-1)*dZ_dXi(p,i-0.5,j-1)*0.25;
                
                if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
                    %dA_dEta = 0.25*(A(i,j) + A(i,j-1) - A(i-1,j) - A(i-1,j-1));
                    S(x,y) = S(x,y) + scale_factor;
                    S(x,y-1) = S(x,y-1) + scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-(p.Xn-1)-1) = S(x,y-(p.Xn-1)-1) - scale_factor;
                    
                elseif(i ==  1) %  top border
                    %dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i,j) - A(i,j-1));
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y-1+(p.Xn-1)) = S(x,y-1+(p.Xn-1)) + scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    
                else % in between top and bottom border. we can use a standard quadrature
                    % dA_dEta = 0.25*(A(i+1,j) + A(i+1,j-1) - A(i-1,j) - A(i-1,j-1));
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                    S(x,y-1+(p.Xn-1)) = S(x,y-1+(p.Xn-1)) + scale_factor;
                    S(x,y-(p.Xn-1)) = S(x,y-(p.Xn-1)) - scale_factor;
                    S(x,y-1-(p.Xn-1)) = S(x,y-1-(p.Xn-1)) - scale_factor;
                end
                
            end
        end
        
        %% Integral 4
        if(true)
            if(i ~=p.Zn-1) % No diffusive flux through the bottom border.
                %dA_dEta = A(i+1,j) - A(i,j);
                %Diffusion
                % output = -2*pi*p.W/(p.Xn-1).*(j-0.5) .* (p.dx * p.dZ_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
                %     dA_dXi_4(i,j,p,A) - 1/( p.dX_dXi_preCalc(i,j,4)*p.dZ_dEta_preCalc(i,j,4))* ...
                %     dA_dEta*(p.dx*p.dZ_dXi_preCalc(i,j,4)^2 + p.dz*p.dX_dXi_preCalc(i,j,4)^2));
                
                
                scale_factor =  -2*pi*p.W/(p.Xn-1).*(j-0.5).*p.dx*dZ_dXi(p,i,j-0.5)/dX_dXi(p,i,j-0.5)*0.25;
                
                if(j == p.Xn-1) % if we are on the right border, we have to change our derivative quadrature accordingly
                    %  dRd_dXi = 0.25*(Rd(i,j) + Rd(i+1,j) - Rd(i,j-1) - Rd(i+1,j-1));
                    
                    S(x,y) = S(x,y) + scale_factor;
                    S(x,y+(p.Xn-1)) =  S(x,y+(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y+(p.Xn-1)-1) = S(x,y+(p.Xn-1)-1) - scale_factor;
                    
                elseif(j ==  1) %  left
                    %  dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i+1,j+1) - Rd(i,j) - Rd(i+1,j));
                    S(x,y+1) = S(x,y+1) + scale_factor;
                    S(x,y+1+(p.Xn-1)) =  S(x,y+1+(p.Xn-1))+ scale_factor;
                    S(x,y) = S(x,y) - scale_factor;
                    S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) - scale_factor;
                    
                else % in between left and right border. we can use a standard quadrature
                    % dRd_dXi = 0.25*(Rd(i,j+1) + Rd(i+1,j+1) - Rd(i,j-1) - Rd(i+1,j-1));
                    S(x,y+1) = S(x,y+1) + scale_factor;
                    S(x,y+1+(p.Xn-1)) =  S(x,y+1+(p.Xn-1))+ scale_factor;
                    S(x,y-1) = S(x,y-1) - scale_factor;
                    S(x,y+(p.Xn-1)-1) = S(x,y+(p.Xn-1)-1) - scale_factor;
                end
                
                scale_factor =  2*pi*p.W/(p.Xn-1).*(j-0.5).*1/( dX_dXi(p,i,j-0.5)*dZ_dEta(p,i,j-0.5))*(p.dx*dZ_dXi(p,i,j-0.5)^2 + p.dz*dX_dXi(p,i,j-0.5)^2);
                S(x,y+(p.Xn-1)) = S(x,y+(p.Xn-1)) + scale_factor;
                S(x,y) = S(x,y) - scale_factor;
                
            end
            
         
        end
        
    end
    
end

% division by alement area
 temp = p.volumes_cyl';
 temp2 = ones(2*(p.Xn-1)*(p.Zn-1) + 2*(p.Xn-1),1);
 temp2(1:length(temp(:))) = temp(:);
 temp2(length(temp(:))+1:2*length(temp(:)) ) = temp(:);
 S = S./temp2;

end

