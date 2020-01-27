
% Internal Nodes
for j = 2:p.Xn-2 % xi
    for i = 2:p.Yn-2 % eta
        
        %%%%%%%% DISSOLVED NUTRIENTS %%%%%%%%%%%
        % Integral 1 - precalc
        dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
            ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i+1,j) + Rd(i+1,j+1) - ...
            Rd(i-1,j) - Rd(i-1,j+1) ));
        
        % Integral 2 - precalc
        dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j+1) + ...
            Rd(i,j+1) - Rd(i-1,j-1) - Rd(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
            (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
        
        % Integral 3 - precalc
        dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
            p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i+1,j-1)+Rd(i+1,j)-Rd(i-1,j-1)-Rd(i-1,j)) - ...
            p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));
        
        % Integral 4 - precalc
        dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
            0.25*(Rd(i+1,j+1) + Rd(i,j+1) - Rd(i,j-1) - Rd(i+1,j-1)) + ...
            1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
            (Rd(i+1,j) - Rd(i,j) )* ...
            (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
        
        test = 1;
    end
end


%%%%%%%%%% BORDERS %%%%%%%%%%%
% We approximate the derivatives on the border using the fact that vhat * (gradient*Rd) = 0,
% where vhat is the transformed normal vector of the face we are
% integrating. Note that vhat is not neccesarily parallell with the axes
% xi,eta as was previously assumed. Using the relation vhat =
% nhat^T * J^-1 * J^-T, where J is the Jacobian matrix, we can calculate
% the components of nhat and use it to approximate derivatives on the
% borders. This results in a centered approximation of all derivatives,
% even those on the border.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%% BOTTOM BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:p.Xn-2 % xi
    i=p.Yn -1; % eta, bottom border
    
    % Integral 1
    dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i,j) + Rd(i,j+1) - ...
        Rd(i-1,j) - Rd(i-1,j+1) ));
    
    % Integral 2 - precalc
    dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j+1) + ...
        Rd(i,j+1) - Rd(i-1,j-1) - Rd(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    
    % Integral 3
    dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i,j-1)+Rd(i,j)-Rd(i-1,j-1)-Rd(i-1,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));
    
    % Integral 4 = 0
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%% LEFT BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:p.Yn-2 % eta
    j=1; % Left border
    % Integral 1
    dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        (p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i+1,j) + Rd(i+1,j+1) - ...
        Rd(i-1,j) - Rd(i-1,j+1) ));
    
    %     %Integral 2
    dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j+1) + ...
        Rd(i,j+1) - Rd(i-1,j) - Rd(i,j)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Integral 3 = 0 due to boundary conditions
    
    % Integral 4 - Rdlternative (correct?) Derivation - precalc
    dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(Rd(i+1,j+1) + Rd(i,j+1) - Rd(i,j) - Rd(i+1,j)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (Rd(i+1,j) - Rd(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
    test = 1;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:p.Yn-2 % eta
    j=p.Xn-1; % xi, right border
    
    % integral 1 = 0 due to bc
    
    % Integral 2
    dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j) + ...
        Rd(i,j) - Rd(i-1,j-1) - Rd(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Integral 3
    dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i+1,j-1)+Rd(i+1,j)-Rd(i-1,j-1)-Rd(i-1,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));
    
    % Integral 4 -  precalc
    dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(Rd(i+1,j) + Rd(i,j) - Rd(i,j-1) - Rd(i+1,j-1)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (Rd(i+1,j) - Rd(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% TOP BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:p.Xn-2 % xi
    i=1; % eta, top border
    
    %     % Integral 1
    dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i+1,j) + Rd(i+1,j+1) - ...
        Rd(i,j) - Rd(i,j+1) ));
    
    % Integral 2 = 0
    
    % % Integral 3
    dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*(p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i+1,j-1)+Rd(i+1,j)-Rd(i,j-1)-Rd(i,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));
    
    % Integral 4 - precalc
    dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(Rd(i+1,j+1) + Rd(i,j+1) - Rd(i,j-1) - Rd(i+1,j-1)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (Rd(i+1,j) - Rd(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
    test= 1;
end


%% %%%%%%%%%%%%%%%% CORNERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%% Top Left %%%% precalcad
i=1; %eta
j=1; %Xi

% Integral 2 = 0
% Integral 3 = 0

% Integral 1
dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i+1,j) + Rd(i+1,j+1) - ...
    Rd(i,j) - Rd(i,j+1) ));

% Integral 4
dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
    0.25*(Rd(i+1,j+1) + Rd(i,j+1) - Rd(i,j) - Rd(i+1,j)) + ...
    1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
    (Rd(i+1,j) - Rd(i,j) )* ...
    (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);

%% %%%% Top Right %%%% precalcad
i=1; %eta
j=p.Xn-1; %Xi

% % Integral 3
dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
    p.dY_dEta_preCalc(i,j,3)) * p.dx*(p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i+1,j-1)+Rd(i+1,j)-Rd(i,j-1)-Rd(i,j)) - ...
    p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));

% Integral 4 - precalc
dRdt(i,j) =  dRdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
    0.25*(Rd(i+1,j) + Rd(i,j) - Rd(i,j-1) - Rd(i+1,j-1)) + ...
    1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
    (Rd(i+1,j) - Rd(i,j) )* ...
    (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);

%% %%%% Bottom Left %%%% precalcad
i=p.Yn-1; %eta
j=1; %Xi

% integral 3 = 0
% integral 4 = 0

% Integral 1
dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i,j) + Rd(i,j+1) - ...
    Rd(i-1,j) - Rd(i-1,j+1) ));

% Integral 1 - boundary value approx
% dRdt(i,j) = dRdt(i,j) +  p.dx * 1/dX_dXi(p,i-0.5,j) * ...
%     ( dY_dEta(p,i-0.5,j) *(Rd(i,j+1) - Rd(i,j)) - dY_dXi(p,i-0.5,j) * ...
%     0.5*( (-1)*abs(1/(dX_dXi(p,i,j)^2*dY_dEta(p,i,j)^2)*dY_dEta(p,i,j)*dY_dXi(p,i,j))*(Rd(i,j+1)-Rd(i,j)) +...
%     0.5*(Rd(i-1,j+1) + Rd(i,j+1) - Rd(i,j) -Rd(i-1,j))) );

% Integral 2
dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j+1) + ...
    Rd(i,j+1) - Rd(i-1,j) - Rd(i,j)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
    (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));

%% %%%% Bottom Right %%%% precalcad
i=p.Yn-1; %eta
j=p.Xn-1; %Xi

% integral 1 = 0
% integral 4 = 0

% % Integral 2
dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j) + ...
    Rd(i,j) - Rd(i-1,j-1) - Rd(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
    (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));

% Integral 3
dRdt(i,j) = dRdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
    p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(Rd(i,j-1)+Rd(i,j)-Rd(i-1,j-1)-Rd(i-1,j)) - ...
    p.dY_dEta_preCalc(i,j,3) * (Rd(i,j) - Rd(i,j-1)));
