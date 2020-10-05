clear
close all

p = struct;

n = 40; % resolution of the grid (>= 2)
p.n = n;
p.diff_coeff = 3;
p.c = 1;
%p.c = 16.66667; % increase results in lower diffusion, given fixed diffusion coeff.

% the origin is in the upper left corner, with the top right corner being
% (1,n), the bottom left (n,1) and the bottom right corner (n,n).
% the flow is oriented from (1,1) to (1,2) with the other flows following suit.

T = 20;
p.T = T; % number of timesteps
p.A0 = ones(p.n);
% p.A0(1,1) = 2;
% p.A0(5,5) = 2;
% p.A0(1,5) = 2;
% p.A0(5,1) = 2;
% p.A0(3,3) = 2;

p.A0((ceil(n/2)-n/10:ceil(n/2)+n/10),(ceil(n/2)-n/10:ceil(n/2)+n/10)) = 2; % setting the concentration in the center to 1 at t = 0
%p.A0(2,2) = 2;
p.v = 2; % speed
p.L = 20; % spacial grid size
% delta_t = 0.0005; 5 length of timestep

% diffusion with diffusion coefficient set to v/L
p.A_diff0 = ones(p.n);
p.A_diff0((ceil(n/2)-n/10:ceil(n/2)+n/10),(ceil(n/2)-n/10:ceil(n/2)+n/10)) = 2; % setting the concentration in the center to 1 at t = 0

%% run simulation
ode_opts = odeset( 'abstol' , 1e-5 , 'reltol' , 1e-5);

%T_vec = [0,ceil(p.T/2), p.T];
T_vec = [0,10, p.T];

[t_circ,A] = ode15s( @(t,Y) rhs_circmodel(t,Y,p), T_vec , p.A_diff0 , ode_opts);
[t_diff,A_diff] = ode15s( @(t,Y) rhs_standard_diff(t,Y,p), T_vec , p.A0 , ode_opts);

test = 1;
%% plot

if(true)
    if(abs(sum(sum(A(3,:)))- sum(sum(A(1,:)))) > 10e-5)
        disp("LEAK!");
    end
    
    if(abs(sum(sum(A_diff(3,:)))- sum(sum(A(1,:)))) > 10e-5)
        disp("LEAK!");
    end
    
    
    
    
    subplot(2,3,1)
    h  = surf(p.A_diff0);
    set(h,'edgecolor','none');
    zlim([1 2])
    
    subplot(2,3,2)
    %t_temp = ceil(length(A_diff(:,1))/2);
    temp = reshape(A_diff(2,:),p.n,p.n); % state at half simulation time
    temp(1,:) = NaN;
    temp(end,:) = NaN;
    temp(:,1) = NaN;
    temp(:,end) = NaN;
  
    h  = surf(temp);
    set(h,'edgecolor','none');
    zlim([1 2])
    
    subplot(2,3,3)
    temp = reshape(A_diff(end,:),p.n,p.n); % end simulation state
    temp(1,:) = NaN;
    temp(end,:) = NaN;
    temp(:,1) = NaN;
    temp(:,end) = NaN;
    h  = surf(temp);
    set(h,'edgecolor','none');
    zlim([1 2])
    
    % ____________________________
    subplot (2,3,4)
    h  = surf(p.A0);
    set(h,'edgecolor','none');
    zlim([1 2])
    
    subplot (2,3,5)
    temp = reshape(A(2,:),p.n,p.n); % state at half simulation time
    temp(1,:) = NaN;
    temp(end,:) = NaN;
    temp(:,1) = NaN;
    temp(:,end) = NaN;
    h  = surf(temp);
    set(h,'edgecolor','none');
    zlim([1 2])
    
    subplot (2,3,6)
    temp = reshape(A(3,:),p.n,p.n); % state at half simulation time
    temp(1,:) = NaN;
    temp(end,:) = NaN;
    temp(:,1) = NaN;
    temp(:,end) = NaN;
    h  = surf(temp);
    set(h,'edgecolor','none');
    zlim([1 2])
end

temp_diff = reshape(A_diff(2,:),p.n,p.n); % state at half simulation time
temp_diff(1,:) = [];
temp_diff(end,:) = [];
temp_diff(:,1) = [];
temp_diff(:,end) = [];

temp_circ = reshape(A(2,:),p.n,p.n); % state at half simulation time
temp_circ(1,:) = [];
temp_circ(end,:) = [];
temp_circ(:,1) = [];
temp_circ(:,end) = [];

deviation = norm((temp_circ -temp_diff)./temp_diff);
disp(deviation);

% deviation = mean((abs(A(2,:)-A_diff(2,:)))./A_diff(2,:));
% disp(deviation);
% deviation = mean((abs(A(3,:)-A_diff(3,:)))./A_diff(3,:));
% disp(deviation);
disp("________________________________");