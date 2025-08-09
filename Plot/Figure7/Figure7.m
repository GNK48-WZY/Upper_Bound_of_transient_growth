% Initialize the cell array
all_eigenvalues = cell(1, 7);
t0 = [0, 10, 20, 40, 60, 80, 100];
dt = 1;
% Loop through files 1 to 7
for k = 1:7
    % Construct filename
    filename = sprintf('all_eigenvalues_%d.csv', k);
    
    % Read data
    cell_data = readcell(filename);
    
    % Convert to numeric matrix
    if iscell(cell_data)
        numeric_vector = cell2mat(cell_data);
    else
        numeric_vector = cell_data;
    end
    
    % Ensure column orientation
    if isrow(numeric_vector)
        numeric_vector = numeric_vector';
    end
    
    % Store numeric matrix
    all_eigenvalues{k} = numeric_vector;
end


%% Plot minimum eigenvalues
figure;
hold on;
cmap = lines(length(t0)); % Colormap
set(groot, 'defaultTextInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex','FontSize',20);

for j3 = 1:length(t0)
    if ~isempty(all_eigenvalues{j3})
        time_points = (1:size(all_eigenvalues{j3},1))*dt;
        % Min across columns for each time point
        min_eig = min(all_eigenvalues{j3}, [], 2);
        plot(time_points, min_eig, 'Color', cmap(j3,:), 'LineWidth', 1.5);
    end
end

% Colorbar
c = colorbar;
colormap(cmap);
c.Ticks = linspace(0, 1, length(t0));
c.TickLabels = arrayfun(@(x) num2str(x), t0, 'UniformOutput', false);
c.FontSize = 20;

legend({'$\lambda_{min}[\mathbf{P}(t)]$'}, 'Location', 'best', ...
    'Interpreter', 'latex','FontSize',24);

xlabel('t','FontSize',20);
grid on; box on;
hold off;

%% Plot maximum eigenvalues
figure;
hold on;
cmap = lines(length(t0));
set(groot, 'defaultTextInterpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize',20);

for j3 = 1:length(t0)
    if ~isempty(all_eigenvalues{j3})
        time_points = (1:size(all_eigenvalues{j3},1))*dt;
        % Max across columns for each time point
        max_eig = max(all_eigenvalues{j3}, [], 2);
        plot(time_points, max_eig, '--', 'Color', cmap(j3,:), 'LineWidth', 1.5);
    end
end

% Colorbar
c = colorbar;
colormap(cmap);
c.Ticks = linspace(0, 1, length(t0));
c.TickLabels = arrayfun(@(x) num2str(x), t0, 'UniformOutput', false);
c.FontSize = 20;

legend({'$\lambda_{max}[\mathbf{P}(t)]$'}, 'Location', 'best', ...
    'Interpreter', 'latex', 'FontSize', 24);

set(gca, 'yscale', 'log');
xlabel('t','FontSize', 20);
grid on; box on;
hold off;
