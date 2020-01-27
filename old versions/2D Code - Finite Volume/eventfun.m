function [x,isterm,dir] = eventfun(t,Y,p,I)


%% reformatting input and calculation of light intensity etc.
% separating the state variables from Y
A = Y(1:(p.Xn-1)*(p.Yn-1));
Rd = Y((p.Xn-1)*(p.Yn-1)+1 : 2*(p.Xn-1)*(p.Yn-1));
Rs = Y(2*(p.Xn-1)*(p.Yn-1) +1 : 2*(p.Xn-1)*(p.Yn-1) + (p.Xn-1));

A = reshape(A,[(p.Xn-1) (p.Yn-1)]);
Rd = reshape(Rd,[(p.Xn-1) (p.Yn-1)]);
A = A';
Rd = Rd';

% Calculation of the light intensity at each grid element center
I = zeros(p.Yn-1,p.Xn-1);
for x = 1:p.Xn-1
    integral = 0;
    for y=1:p.Yn-1
        %y_coord = 0.25*(p.Y(y,x) + p.Y(y+1,x) + p.Y(y,x+1) + p.Y(y+1,x+1));
        integral = integral + p.Y_vol(y,x) *A(y,x);
        I(y,x) = p.I0 * exp(-p.k*integral -p.kbg*p.Y_vol(y,x));
    end
end

% Calculation of the algal production G(R,I)
G =  p.Gmax .* Rd./(Rd+p.M) .* I./(I+p.H);

%% Assembling the right hand side of our eq system
dAdt = zeros(p.Yn-1, p.Xn-1);

%% diffusion & convection

% Using the divergence theorem, the surface integral of a square control volume in
% our mesh is expressed as the line integral along the border of the control
% volume. This line integral is split up into four parts, computing the integral
% along each side separately. The integrals are numbered in a counter-clockwise fashion,
% where integral 1 is the integral from the bottom right corner to the top right corner.

div_area = true;

% Internal Nodes
for j = 2:p.Xn-2 % xi
    for i = 2:p.Yn-2 % eta
        
        %%%%%%%%% ALGAE - DIFFUSION %%%%%%%%%%%%%
        % Integral 1 - precalc
        dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
            ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i+1,j) + A(i+1,j+1) - ...
            A(i-1,j) - A(i-1,j+1) ));
        
        % Integral 2 - precalc
        dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j+1) + ...
            A(i,j+1) - A(i-1,j-1) - A(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
            (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
        
        % Integral 3 - precalc
        dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
            p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i-1,j-1)-A(i-1,j)) - ...
            p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));
        
        % Integral 4 - precalc
        dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
            0.25*(A(i+1,j+1) + A(i,j+1) - A(i,j-1) - A(i+1,j-1)) + ...
            1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
            (A(i+1,j) - A(i,j) )* ...
            (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
        
                        %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));
        
    end
end


%%%%%%%%%% BORDERS %%%%%%%%%%%
% We approximate the derivatives on the border using the fact that vhat * (gradient*A) = 0,
% where vhat is the transformed normal vector of the face we are
% integrating. Note that vhat is not neccesarily parallell with the axes
% xi,eta as was previously assumed. Using the relation vhat =
% nhat^T * J^-1 * J^-T, where J is the Jacobian matrix, we can calculate
% the components of nhat and use it to approximate derivatives on the
% borders. This results in a centered approximation of all derivatives,
% even those on the border.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:p.Xn-2 % xi
    i=p.Yn -1; % eta, bottom border
    
    %%%%%%%%%%%%% ALGAE - DIFFUSION %%%%%%%%%%%%%%%
    % Integral 1
    dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i,j) + A(i,j+1) - ...
        A(i-1,j) - A(i-1,j+1) ));

    % Integral 2 - precalc
    dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j+1) + ...
        A(i,j+1) - A(i-1,j-1) - A(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Integral 3
    dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(A(i,j-1)+A(i,j)-A(i-1,j-1)-A(i-1,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));

    
    % Integral 4 = 0
    
    
    %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Algae sink with a 
    
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) =  dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
    
      %Integral 4 = 0 on the bottom
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% LEFT BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:p.Yn-2 % eta
    j=1; % Left border
    
    %%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
    % Integral 1
    dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        (p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i+1,j) + A(i+1,j+1) - ...
        A(i-1,j) - A(i-1,j+1) ));
    
    %     %Integral 2
    dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j+1) + ...
        A(i,j+1) - A(i-1,j) - A(i,j)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Integral 3 = 0 due to boundary conditions
    
    % Integral 4 - Alternative (correct?) Derivation - precalc
    dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(A(i+1,j+1) + A(i,j+1) - A(i,j) - A(i+1,j)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (A(i+1,j) - A(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
        %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) -  p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% RIGHT BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:p.Yn-2 % eta
    j=p.Xn-1; % xi, right border
    
    %%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
    % integral 1 = 0 due to bc
    % Integral 2
    dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j) + ...
        A(i,j) - A(i-1,j-1) - A(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    % Integral 3
    dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i-1,j-1)-A(i-1,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));
    
    %     % Integral 4
    %     dAdt(i,j) = dAdt(i,j) + 1/abs(dX_dXi(p,i,j-0.5) * dY_dEta(p,i,j-0.5))*p.dy*dX_dXi(p,i,j-0.5)^2 * ...
    %         (A(i+1,j)-A(i,j)) + p.dx*dY_dXi(p,i,j-0.5)^2*(A(i+1,j)-A(i,j)) -  ...
    %         p.dx*dY_dXi(p,i,j-0.5)*dY_dEta(p,i,j-0.5)*0.25*(A(i,j) + A(i+1,j) - A(i,j-1) - A(i+1,j-1));
    %
    % Integral 4 -  precalc
    dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(A(i+1,j) + A(i,j) - A(i,j-1) - A(i+1,j-1)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (A(i+1,j) - A(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
            %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) -  p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% TOP BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:p.Xn-2 % xi
    i=1; % eta, top border
    
    %%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
    %     % Integral 1
    dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i+1,j) + A(i+1,j+1) - ...
        A(i,j) - A(i,j+1) ));

    % Integral 2 = 0
    
    % % Integral 3
    dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*(p.dY_dXi_preCalc(i,j,3)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i,j-1)-A(i,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));
    
    % Integral 4 - precalc
    dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
        0.25*(A(i+1,j+1) + A(i,j+1) - A(i,j-1) - A(i+1,j-1)) + ...
        1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
        (A(i+1,j) - A(i,j) )* ...
        (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);
    
                %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2 = 0  from bc
    %dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));

    test= 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CORNERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Top Left %%%% precalcad
i=1; %eta
j=1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% Integral 2 = 0
% Integral 3 = 0

% Integral 1
dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i+1,j) + A(i+1,j+1) - ...
    A(i,j) - A(i,j+1) ));

% Integral 4
dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
    0.25*(A(i+1,j+1) + A(i,j+1) - A(i,j) - A(i+1,j)) + ...
    1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
    (A(i+1,j) - A(i,j) )* ...
    (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);

                %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2 = 0 from bc
    %dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));

%% %%%% Top Right %%%% precalcad
i=1; %eta
j=p.Xn-1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% % Integral 3
dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
    p.dY_dEta_preCalc(i,j,3)) * p.dx*(p.dY_dXi_preCalc(i,j,3)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i,j-1)-A(i,j)) - ...
    p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));

% Integral 4 - precalc
dAdt(i,j) =  dAdt(i,j) - p.dx * p.dY_dXi_preCalc(i,j,4)/p.dX_dXi_preCalc(i,j,4)* ...
    0.25*(A(i+1,j) + A(i,j) - A(i,j-1) - A(i+1,j-1)) + ...
    1/( p.dX_dXi_preCalc(i,j,4) *p.dY_dEta_preCalc(i,j,4))*...
    (A(i+1,j) - A(i,j) )* ...
    (p.dx*p.dY_dXi_preCalc(i,j,4)^2 + p.dy*p.dX_dXi_preCalc(i,j,4)^2);

                %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2 = 0 from bc
    %dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*0.5*(A(i,j) + A(i+1,j));

%% %%%% Bottom Left %%%% precalcad
i=p.Yn-1; %eta
j=1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% integral 3 = 0
% integral 4 = 0

% Integral 1
dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i,j) + A(i,j+1) - ...
    A(i-1,j) - A(i-1,j+1) ));

% Integral 2
dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j+1) + ...
    A(i,j+1) - A(i-1,j) - A(i,j)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
    (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));

    %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Algae sink with a 
    
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4 = 0 on the bottom
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);

%% %%%% Bottom Right %%%% precalcad
i=p.Yn-1; %eta
j=p.Xn-1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% integral 1 = 0
% integral 4 = 0

% % Integral 2
dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j) + ...
    A(i,j) - A(i-1,j-1) - A(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
    (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));

% Integral 3
dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
    p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(A(i,j-1)+A(i,j)-A(i-1,j-1)-A(i-1,j)) - ...
    p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));

    %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Algae sink with a 
    
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) = dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
      %Integral 4 = 0 on the bottom
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);

%% Divion by element area
% We want the rate of change of the concentration, but the above
% calculations are in total amount of algae/nutrients. We divide by the are
% of each element to yield concentrations.

if(div_area)
    dAdt = dAdt./p.vol_areas;
end

%% algal net growth
dAdt = dAdt + (G - p.lbg).*A;

    
   %% Evaluating size of derivatives
x = norm(dAdt) - 1e-5;
isterm = 1;
dir = -1;  %or -1, doesn't matter
end