clear all;

all_V_lyap = readmatrix("all_V_lyap.csv");
bound = readmatrix("bound.csv");
time_points = readmatrix("time_points.csv");
t_sim = readmatrix("t_sim.csv");
cmap = lines(5);
figure;
hold on;

% Plot all three trajectories with the same legend entry
h_traj_group = gobjects(1, 3);  % Handle array for trajectory group
for traj_idx = 1:3
    h = semilogy(t_sim, all_V_lyap(traj_idx, :), '--', 'LineWidth', 1.5, ...
                 'Color', cmap(traj_idx, :), ...
                 'HandleVisibility', 'off');  % Hide individual legend entries
    
    h_traj_group(traj_idx) = h;
end

% Create a single proxy line for legend (using first trajectory's properties)
h_proxy = semilogy(NaN, NaN, '--', 'LineWidth', 1.5, 'Color', 'black', ...
                  'DisplayName', '$\mathbf{x}^{\ast}(t)\mathbf{P}(t)\mathbf{x}(t)$');

% Plot bound
h_bound = semilogy(time_points, bound, '-', 'LineWidth', 2, 'Color', cmap(4,:), ...
                 'DisplayName', '$\lambda_{\max}[\mathbf{P}(t)]G(t)$');

xlim([0 200]);
set(gca, 'YScale', 'log');

% Configure legend with proxy and bound
legend([h_proxy, h_bound], 'Interpreter', 'latex', 'Location', 'best');

% Axis labels and formatting
xlabel('$t$', 'Interpreter', 'latex', 'FontSize',16);
title('');
ylabel('');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
box on;
grid on;

