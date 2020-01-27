function [Ma,Mrb,M] = Stiffness_matrix(p)

%% Mass matrix of the linearized system
% assembles the mass matrix of the linearized system with integrated
% boundary conditions.


M = zeros(2*(p.Xn-1)*(p.Zn-1) + 2*(p.Xn-1));


%% mass matrix for the algae
Ma = zeros((p.Xn-1)*(p.Zn-1));

for j = 1:p.Xn-1 % xi
    for i = 1:p.Zn-1 % eta
        % Diffusion and sinking of algae with  BC:s
        dAdt(i,j) = integral_A_1(i,j,p,A) + integral_A_2(i,j,p,A) + ...
            integral_A_3(i,j,p,A) + integral_A_4(i,j,p,A);
        
        % Integral 1
        if(j ~= p.Xn-1) % we are not on the right border.
            
            %%%%%%%%% remove when done, only kept for reference %%%%%%%%%%
            M(i,j) = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
                (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi - ...
                p.4dZ_dXi_preCalc(i,j,1) * dA_dEta_1(i,j,p,A));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            444444
            
            %  implementation of 2*pi*p.W/(p.Xn-1).*j .* p.dx * 1/p.dX_dXi_preCalc(i,j,1) * (p.dZ_dEta_preCalc(i,j,1) .* dA_dXi
            % into the mass matrix.
            % A simple forward difference is used to approximate dA_dXi,
            % dA_dXi =  A(i,j+1) - A(i,j).
            
            M(i,j) = - 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi(i,j,1) * ...
                p.dZ_dEta(i,j,1);
            
            M(i,j+1) = 2*pi*p.W/(p.Xn-1).*j .*  p.dx * 1/p.dX_dXi(i,j,1) * ...
                p.dZ_dEta(i,j,1);
            
            
            if(i ==  p.Zn-1) % if we are on the bottom border, we have to change our derivative quadrature accordingly
                dA_dEta = 0.25*(A(i,j) + A(i,j+1) - A(i-1,j) - A(i-1,j+1));                
            elseif(i ==  1) %  top border
                dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i,j) - A(i,j+1));
            else % in between top and bottom border. we can use a standard quadrature
                dA_dEta = 0.25*(A(i+1,j) + A(i+1,j+1) - A(i-1,j) - A(i-1,j+1));
            end
            
            
        end
        
        
        % BC:s imply that there is no flux on the right border, and the mass
        % Matrix values corresponding to those nodes is simply zero.
        
        
    end
    
end




%% mass matrix for the dissolved nutrients
Mrb = zeros((p.Xn-1)*(p.Zn-1));

end

