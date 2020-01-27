function [output] = dAdt_fun(t,Y,p,I)
% Returns the right hand side of the equation dAdt = r.h.s,
% where the output is a column vector and the nodes are ordered
% row by row, so that the first p.Xn nodes are the A values of the
% surface layer and so forth.
%
% function is used as input to ode15s

% The vectors I,A,R,Rs are allocated in the main script and passed as input
% to avoid allocation each time the function is called.

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
dRdt = zeros(p.Yn-1, p.Xn-1);
dRsdt = zeros(1,p.Xn-1);

%% diffusion & convection

% Using the divergence theorem, the surface integral of a square control volume in
% our mesh is expressed as the line integral along the border of the control
% volume. This (closed) line integral is split up into four parts, computing the integral
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
        
        %%%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%
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
        
        %test = 1;
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
    
    % Integral 1 - boundary value approximation
    %             dAdt(i,j) = dAdt(i,j) -  p.dx * 1/dX_dXi(p,i-0.5,j) * ...
    %         ( dY_dEta(p,i-0.5,j) *(A(i,j+1) - A(i,j)) - dY_dXi(p,i-0.5,j) * ...
    %         0.5*( (-1)*abs(1/(dX_dXi(p,i,j)^2*dY_dEta(p,i,j)^2)*(dY_dXi(p,i,j)^2+dX_dXi(p,i,j)^2  ))*(A(i,j+1)-A(i,j)) +...
    %         0.5*(A(i-1,j+1) + A(i,j+1) - A(i,j) -A(i-1,j))) );
    
    % Integral 2 - precalc
    dAdt(i,j) = dAdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(A(i-1,j+1) + ...
        A(i,j+1) - A(i-1,j-1) - A(i,j-1)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
        (A(i,j)-A(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));
    
    
    % Integral 3
    dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
        p.dY_dEta_preCalc(i,j,3)) * p.dx*( p.dY_dXi_preCalc(i,j,3)*0.25*(A(i,j-1)+A(i,j)-A(i-1,j-1)-A(i-1,j)) - ...
        p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));
    
    % Integral 3 - boundary value approximation
    %     dAdt(i,j) = dAdt(i,j) + sqrt( dX_dEta(p,i-0.5,j-1)^2 + dY_dEta(p,i-0.5,j-1)^2) / abs( dX_dXi(p,i-0.5,j-1)* ...
    %         dY_dEta(p,i-0.5,j-1)) * p.dx*( dY_dXi(p,i-0.5,j-1)*0.25*(A(i,j-1)+A(i,j)-A(i-1,j-1)-A(i-1,j)) - ...
    %         dY_dEta(p,i-0.5,j-1) * (A(i,j) - A(i,j-1)));
    
    % Integral 4 = 0
    
    
    %%%%%%%%%%%%% ALGAE - ADVECTION %%%%%%%%%%%%%%%
    % Algae sink with a 
    
    % Integral 1 = 0
    % Integral 3 = 0
    %Integral 2
    dAdt(i,j) =  dAdt(i,j) + p.v*p.dX_dXi_preCalc(i,j,2)*0.5*(A(i,j) + A(i-1,j));
    
      %Integral 4 = 0 on the bottom
    dAdt(i,j) = dAdt(i,j) - p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);
    dRsdt(j) = dRsdt(j) + p.q*p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);
    
    %%%%%%%%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%%%%
    
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
    
    
    % Integral 2
    %     dAdt(i,j) = dAdt(i,j) + (p.dx * dY_dXi(p,i-1,j-0.5)/dX_dXi(p,i-1,j-0.5)* ...
    %         0.5*(0.5*(A(i-1,j+1) + A(i,j+1) - A(i-1,j) - A(i,j)) + (A(i,j) -A(i-1,j))* ...
    %         (-1)*abs((dY_dXi(p,i-1,j-1)^2 + dX_dXi(p,i-1,j-1)^2)/(dX_dXi(p,i-1,j-1)^2 * dY_dEta(p,i-1,j-1)^2))) + ...
    %         (-1)*( dX_dXi(p,i-1,j-0.5)*dY_dEta(p,i-1,j-0.5))* ...
    %         (A(i,j)-A(i-1,j))*(p.dx*dY_dXi(p,i-1,j-0.5)^2 + p.dy*dX_dXi(p,i-1,j-0.5)^2));
    
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
    
    %%%%%%%%%%%%%%  DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%%%
    
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
    
    
    %test = 1;
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
    
    %%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%%%
    
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% TOP BORDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 2:p.Xn-2 % xi
    i=1; % eta, top border
    
    %%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
    %     % Integral 1
    dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
        ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i+1,j) + A(i+1,j+1) - ...
        A(i,j) - A(i,j+1) ));
    
    %     dAdt(i,j) = dAdt(i,j) +  p.dx * 1/dX_dXi(p,i-0.5,j) * ...
    %         ( dY_dEta(p,i-0.5,j) *(A(i,j+1) - A(i,j)) - dY_dXi(p,i-0.5,j) * ...
    %         0.5*( (A(i,j+1)-A(i,j)) *(-1)*abs(1/(dX_dXi(p,i-1,j)^2*dY_dEta(p,i-1,j)^2)*dY_dEta(p,i-1,j)*dY_dXi(p,i-1,j))+...
    %         0.5*( A(i+1,j)+A(i+1,j+1) - A(i,j) - A(i,j+1) ) )     );
    
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
    
    %%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%
    
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
    
    %test= 1;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CORNERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Top Left Corner %%%% precalcad
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
% dAdt(i,j) = dAdt(i,j) + 1/abs(dX_dXi(p,i,j-0.5) * dY_dEta(p,i,j-0.5))*p.dy*dX_dXi(p,i,j-0.5)^2 * ...
%     (A(i+1,j)-A(i,j)) + p.dx*dY_dXi(p,i,j-0.5)^2*(A(i+1,j)-A(i,j)) -  ...
%     p.dx*dY_dXi(p,i,j-0.5)*dY_dEta(p,i,j-0.5)*0.25*(A(i,j+1) + A(i+1,j+1) - A(i,j) - A(i+1,j));

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


%%%%%%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%

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

%% %%%% Top Right Corner %%%% precalcad
i=1; %eta
j=p.Xn-1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% % Integral 3
dAdt(i,j) = dAdt(i,j) + sqrt( p.dX_dEta_preCalc(i,j,3)^2 + p.dY_dEta_preCalc(i,j,3)^2) / abs( p.dX_dXi_preCalc(i,j,3)* ...
    p.dY_dEta_preCalc(i,j,3)) * p.dx*(p.dY_dXi_preCalc(i,j,3)*0.25*(A(i+1,j-1)+A(i+1,j)-A(i,j-1)-A(i,j)) - ...
    p.dY_dEta_preCalc(i,j,3) * (A(i,j) - A(i,j-1)));

% Integral 4
% dAdt(i,j) = dAdt(i,j) + 1/abs(dX_dXi(p,i,j-0.5) * dY_dEta(p,i,j-0.5))*p.dy*dX_dXi(p,i,j-0.5)^2 * ...
%     (A(i+1,j)-A(i,j)) + p.dx*dY_dXi(p,i,j-0.5)^2*(A(i+1,j)-A(i,j)) -  ...
%     p.dx*dY_dXi(p,i,j-0.5)*dY_dEta(p,i,j-0.5)*0.25*(A(i,j) + A(i+1,j) - A(i,j-1) - A(i+1,j-1));


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


%%%%%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%

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

%% %%%% Bottom Left Corner %%%% precalcad
i=p.Yn-1; %eta
j=1; %Xi
%%%%%%%%%%  ALGAE - DIFFUSION %%%%%%%%%%%%%%
% integral 3 = 0
% integral 4 = 0

% Integral 1
dAdt(i,j) = dAdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(A(i,j+1) - A(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( A(i,j) + A(i,j+1) - ...
    A(i-1,j) - A(i-1,j+1) ));

% Integral 1 - boundary value approx
% dAdt(i,j) = dAdt(i,j) +  p.dx * 1/dX_dXi(p,i-0.5,j) * ...
%     ( dY_dEta(p,i-0.5,j) *(A(i,j+1) - A(i,j)) - dY_dXi(p,i-0.5,j) * ...
%     0.5*( (-1)*abs(1/(dX_dXi(p,i,j)^2*dY_dEta(p,i,j)^2)*dY_dEta(p,i,j)*dY_dXi(p,i,j))*(A(i,j+1)-A(i,j)) +...
%     0.5*(A(i-1,j+1) + A(i,j+1) - A(i,j) -A(i-1,j))) );

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
    dRsdt(j) = dRsdt(j) + p.q*p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);

%%%%%%%%%%  DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%%

% integral 3 = 0
% integral 4 = 0

% Integral 1
dRdt(i,j) = dRdt(i,j) +  p.dx * 1/p.dX_dXi_preCalc(i,j,1) * ...
    ( p.dY_dEta_preCalc(i,j,1) *(Rd(i,j+1) - Rd(i,j)) - p.dY_dXi_preCalc(i,j,1) * 0.25 *( Rd(i,j) + Rd(i,j+1) - ...
    Rd(i-1,j) - Rd(i-1,j+1) ));

% Integral 2
dRdt(i,j) = dRdt(i,j) + (p.dx * p.dY_dXi_preCalc(i,j,2)/p.dX_dXi_preCalc(i,j,2)*0.25*(Rd(i-1,j+1) + ...
    Rd(i,j+1) - Rd(i-1,j) - Rd(i,j)) - 1/( p.dX_dXi_preCalc(i,j,2)*p.dY_dEta_preCalc(i,j,2))* ...
    (Rd(i,j)-Rd(i-1,j))*(p.dx*p.dY_dXi_preCalc(i,j,2)^2 + p.dy*p.dX_dXi_preCalc(i,j,2)^2));

%% %%%% Bottom Right Corner %%%% precalcad
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
    dRsdt(j) = dRsdt(j) + p.q*p.v*p.dX_dXi_preCalc(i,j,4)*A(i,j);

%%%%%%%%%%%%% DISSOLVED NUTRIENTS - DIFFUSION %%%%%%%%%%%%


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

%% Divion by element area
% We want the rate of change of the concentration, but the above
% calculations are in total amount of algae/nutrients per grid cell.
% We divide by the area of each element to yield concentrations.

if(div_area)
    dAdt = dAdt./p.vol_areas;
    dRdt = dRdt./p.vol_areas;
    dRsdt = dRsdt./p.L_bottom';
end

%% algal net growth
dAdt = dAdt + (G - p.lbg).*A;

%% Dissolved nutrients
% Net change in dissolved nutrient
dRdt = dRdt - p.q.*(G-p.lbg).*A;

% Here the boundary condition is wrong. The RHS should be divided by the
% diffusion coefficient in the n-hat direction, not dy.

% OLD (wrong)
%dRdt(end,:) = dRdt(end,:) + (p.r/p.dy).*Rs(:)'.*(p.L_bottom'./p.vol_areas(end,:)) ;

% New:
angle = atan((p.Lmax-p.Lmin)/p.W);
diffusion = sqrt((sin(angle)*p.dx)^2 + (cos(angle)*p.dy)^2); % diffusion coefficient pointing in n-hat dir
dRdt(end,:) = dRdt(end,:) + (p.r/(diffusion).*Rs(:)'.*(p.L_bottom'./p.vol_areas(end,:)));

%% Sediment layer
% Sediment remineralizes into the lake as dissolved nutrients

%OLD (wrong!)
%dRsdt = dRsdt -(p.r/p.dy).*Rs';

% Erroneous implementation of the boundary condition!
dRsdt = dRsdt -(p.r/diffusion).*Rs';

%% reshaping matrices for output
% transpose of matrices in order to use colon notation to reshape to vector form.
dAdt = dAdt';
dRdt = dRdt';


output = [dAdt(:); dRdt(:); dRsdt(:)];

end

