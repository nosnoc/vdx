import casadi.*
close all;
clear all;

prob = vdx.Problem();

n_masses = 1000;

% bounds
lbx = -inf; ubx = inf;
lby = -inf; uby = inf;

% constants
k = 10000;   % Spring constant
g = 9.8; % Gravity
d = 15/n_masses;
m = 20/n_masses;   % mass
for ii=1:n_masses
    % x,y coords of each point
    x = SX.sym(['x_' num2str(ii)], 1);
    y = SX.sym(['y_' num2str(ii)], 1);
    pos = [x;y];
    init_x = -5 + ii/n_masses*10;
    init_y = 0;

    % Add position to decision vector
    prob.w.pos(ii) = {pos, [lbx;lby], [ubx;uby], [init_x;init_y]}; % ith mass position
end

% additional ground constraints
for ii=1:n_masses
    pos = prob.w.pos(ii);
    x = pos(1); y = pos(2);
    
    prob.g.hill(ii) = {y+((x-1)^2)-3, 0, inf};
    prob.g.ground(ii) = {y, 0, inf};
    prob.g.slant(ii) = {y+x + 1, 0, inf};
end

% minimize energy
pos = prob.w.pos(1);
prob.f = m*g*pos(2);
for ii=2:n_masses
    pos = prob.w.pos(ii);
    x = pos(1); y = pos(2);
    prev_pos = prob.w.pos(ii-1);
    hyp_prev = norm(pos-prev_pos);

    % Energy = 0.5*kd^2 + mgh
    prob.f = prob.f + 0.5*k*(hyp_prev-d).^2 + m*g*y;;
end

% fix first and last positions
p1 = [-5;5];
prob.w.pos(1).init = p1;
prob.w.pos(1).lb = p1;
prob.w.pos(1).ub = p1;

p_end = [5;5];
prob.w.pos(end).init = p_end;
prob.w.pos(end).lb = p_end;
prob.w.pos(end).ub = p_end;

casadi_opts = struct;

prob.create_solver(casadi_opts);

prob.solve();

%% Plots
% Plot chain
figure
pos = prob.w.pos.res;
hold on
plot(pos(1,:), pos(2,:), 'b-o')
fplot(@(x) 3-(x-1).^2, [-6, 6], 'r--')
fplot(@(x) -x - 1, [-6, 6], 'r--')
yline(0, 'r--')
xlim([-6,6]);
ylim([-1, 6])
hold off

% plot forces applied by the ground
figure
slant_multipliers = prob.g.slant.mult;
hill_multipliers = prob.g.hill.mult;
ground_multipliers = prob.g.ground.mult;
subplot(3,1,1) 
plot(pos(1,:), -ground_multipliers)
subplot(3,1,2)
plot(pos(1,:), -hill_multipliers)
subplot(3,1,3)
plot(pos(1,:), -slant_multipliers)

% force applied to the ground
figure 
plot(pos(1,:), ground_multipliers + hill_multipliers + slant_multipliers)
