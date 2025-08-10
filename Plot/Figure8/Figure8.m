% Remember to unzip all_eigenvectors_1 to 3 
V = readmatrix("V.csv");
t0 = [0, 10, 20, 40, 60, 80, 100];
dt = 1;
N = 32;
a = 1.2; 
b = 0;
[x, D, D2, D4,w] =finitediff(N,2);
expected_len = 2*(N-2);   % 60 for N=32
raw = all_eigenvectors{j3}{j1_max};

% unwrap nested 1x1 cells if present
while iscell(raw) && numel(raw) == 1
    raw = raw{1};
end

% if textual, try converting to numeric
if ischar(raw) || isstring(raw)
    raw = str2num(raw); %#ok<ST2NM>
end

if ~isnumeric(raw)
    error('Unexpected format in all_eigenvectors{%d}{%d}: %s', j3, j1_max, class(raw));
end

% pick the correct column/row/vector as coefficients
if isvector(raw) && length(raw) == expected_len
    coeffs = raw(:);
elseif size(raw,1) == expected_len && size(raw,2) >= max_idx
    coeffs = raw(:, max_idx);        % common case: modes are columns
elseif size(raw,2) == expected_len && size(raw,1) >= max_idx
    coeffs = raw(max_idx, :)';       % modes stored as rows
else
    error('Cannot find a length-%d vector in raw (size = [%d %d]).', expected_len, size(raw,1), size(raw,2));
end

% Map coefficients to physical state q (prefer V* or V\ over inv(V))
if length(coeffs) == expected_len
    q = coeffs;                      % already physical-state vector
elseif isequal(size(V), [expected_len expected_len])
    % try V * coeffs (modal->physical) first, then V\coeffs
    q_try = [];
    if size(V,2) == length(coeffs)
        q_try = V * coeffs;
    end
    if isempty(q_try) && size(V,1) == length(coeffs)
        q_try = V \ coeffs;
    end
    if isempty(q_try) || ~isvector(q_try) || length(q_try) ~= expected_len
        error('Unable to map coeffs (len=%d) to expected_len=%d using V.', length(coeffs), expected_len);
    end
    q = q_try;
else
    error('V has unexpected size [%d %d].', size(V,1), size(V,2));
end

q = q(:);   % ensure column vector
% ----- end replacement block -----



for j3 = 1:length(t0)
    if ~isempty(all_eigenvectors{j3})
               % Full spatial grid (including boundaries)
        Nx = 32;
        kx = a; kz = b;
        k2 = kx^2 + kz^2;
        
        
       % After getting eigenvector q
       for j1_ind = 1:length(t0)
        max_eig(j1_ind) = max(all_eigenvalues{j3}{j1_ind});
       end
        j1_max = 1;
        [~,max_idx] = max(all_eigenvalues{j3}{j1_max});
        t_current = t0(j3) + (j1_max-1)*dt;

        
        
        q = all_eigenvectors{j3}{j1_max}(:, max_idx);
        q = inv(V)*q;
        % Split state vector
        v_interior = q(1:(N-2));
        wy_interior = q((N-2)+1:end);
        
        % Pad with boundary zeros (now size N)
        v_full = [0; v_interior; 0];
        wy_full = [0; wy_interior; 0];
        
        % Compute ∂v/∂y at INTERIOR points
        Dv_interior = D * v_interior;  % D: (N-2)x(N-2) matrix
        Dv_full = [0; Dv_interior; 0]; % Pad boundaries
        
        % Compute velocity components (now all size N)
        u_hat = (1/k2) * (1i*kx * Dv_full - 1i*kz * wy_full);
        v_hat = v_full;
        w_hat = (1/k2) * (1i*kz * Dv_full + 1i*kx * wy_full);
        
        y_grid = linspace(-1, 1, N)';
        x_grid = linspace(0, 2*pi/a, Nx)';
        [X, Y] = meshgrid(x_grid, y_grid);
        
        % Reconstruct physical space perturbation
        u_physical = real(u_hat .* exp(1i*a*X));
        v_physical = real(v_hat .* exp(1i*a*X));
        w_physical = imag(w_hat .* exp(1i*a*X));
        max_abs = max(abs(u_physical(:)));
%         c_limits = [-max_abs, max_abs];
%         c_limits = [-3,3];
        U_interior = U_yi_function(t_current);
        g_current = exp(-k0 * t_current);
        U_full   = [-g_current; U_interior(:); +g_current];  % size N×1
        x_pos    = pi/a;                                    % center in x
        u_shift  = U_full + x_pos;                          % shifted x‐coordinate
        
        figure(j3); clf;
        set(gcf,'Position',[100 100 800 600]);
        % 1) Primary axes for the contour:
        ax1 = axes('Position',[0.10 0.10 0.75 0.8]);
        contourf(ax1, X, Y, u_physical, 20, 'LineColor','none');
        colormap(ax1, bluewhitered);
        colorbar(ax1,'Location','eastoutside');
        caxis(ax1,[-2 2]);
        set(ax1,'FontSize',26,'Box','on');
        
        



        % 2) Overlay a second axes for U(y):
        ax2 = axes('Position',get(ax1,'Position'), ...
                   'XAxisLocation','top', ...
                   'YAxisLocation','right', ...
                   'Color','none', ...
                   'YTick',[], ...      % hide its y‑ticks
                   'Box','off', ...
                   'XColor','#77AC30', ...
                   'FontSize',26);
        set(ax2, 'YColor', 'none');      % <— hide that vertical line

        
        % 3) Link Y‑limits so the laminar curve sits at the correct vertical positions:
        linkprop([ax1 ax2],{'YLim'});
        ylim(ax2, ax1.YLim);
        
        % 4) Give ax2 its own X‑limits = [-1,1] (laminar velocity range):
        xlim(ax2, [-1 1]);
        
        % 5) Force equal data‑unit aspect on ax2:
        %    one “unit” in x (e.g. from -1 to 0) equals one “unit” in y
        daspect(ax2, [1 1 1]);
        
        % 6) Set evenly‑spaced ticks from -1 to +1:
        lam_ticks = linspace(-1,1,5);
        set(ax2, 'XTick', lam_ticks, 'XTickLabel', arrayfun(@num2str, lam_ticks,'Uni',0));
        % after creating ax2 …
        set(ax2, 'Units', 'centimeters');          % work in cm
        pos = get(ax2, 'Position');                % [left bottom width height]
        pos(1) = pos(1) - 0.75;                       % move left by 1 cm
        set(ax2, 'Position', pos, 'Units', 'normalized');  % apply and switch back
%         axis(ax2,'equal');            

        
        
        % 7) Plot the laminar profile on ax2:
        hold(ax2,'on');
        plot(ax2, U_full, y_grid, 'LineWidth',3, 'Color','#77AC30');
        hold(ax2,'off');

    end
end



function [x, D1, D2, D4,w] =finitediff(N,L)

%Differentiation matrix using finite difference scheme.
%This is suitable for Dirichlet boundary condition v(x=L/2)=v(x=-L/2)=0 at
%the boundary and Neumann boundary condition v'(x=L/2)=v'(x=-L/2)=0. 

dx=L/N;%get grid spacing. 

x=linspace(-L/2,L/2,N);%get the grid point location

x=x(2:end-1); %delete the first and last points that are zero. 

w=dx*ones(1,N-2);%integration weighting 

N_diff=N-2;%The size of differentiation matrices

%First order derivative based on central difference
%f'(x_i)=(f(x_{i+1})-f(x_{i-1})/(2*dx)
%We also use the boundary condition that x_0=x_N=0
D1_stencil=diag(-1*ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1);
D1=D1_stencil/(2*dx);

%Second order derivative based on central difference
%f''(x_i)=(f(x_{i+1})-2f(x_i)+f(x_{i-1})/(dx^2)
%We also use the boundary condition that x_0=x_N=0
D2_stencil=(diag(-2*ones(1,N_diff))+diag(ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1));
D2=D2_stencil/dx^2;

%Forth order derivative based on central difference
%f''''(x_i)=(f(x_{i+2})-4f(x_{i+1})+6f(x_i)-4f(x_{i-1})+f(x_{i-2})/(dx^4)
%This differentiation matrix only go through x_1 up to x_{N-1}
D4_stencil=(diag(6*ones(1,N_diff))+diag(-4*ones(1,N_diff-1),1)+diag(-4*ones(1,N_diff-1),-1)...
    + diag(ones(1,N_diff-2),2)+diag(ones(1,N_diff-2),-2));

%Here, we use the Neumann boundary condition that x_0'=x_N'=0
%such that x_1=x_{-1} and x_{N-1}=x_{N+1} for the ghost points. Then we
%also use the condition x_0 and x_N=0 to express all values based on x_1 up
%to x_{N-1}
D4_stencil(1,1)=D4_stencil(1,1)+1;
D4_stencil(end,end)=D4_stencil(end,end)+1;

D4=D4_stencil/dx^4;

end
