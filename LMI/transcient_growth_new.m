%acc(a, b, k, Re)
% cluster=1;
% if cluster
%     %n_samples=12;
%    % parpool(12);
%    addpath(genpath('/home/zhw22003/YALMIP-master'))
%    system('export PATH=$PATH:/home/zhw22003/mosek/10.2/tools/platform/linux64x86/bin')
%    addpath(genpath('/home/zhw22003/mosek'))
% %     addpath(genpath('/home/jtw24004/YALMIP-master'));
% %     system('export PATH=$PATH:/home/jtw24004/mosek/10.2/tools/platform/linux64x86/bin')
% %     addpath(genpath('/home/jtw24004/mosek'));
% %     addpath(genpath('/home/jtw24004/matlabfunction'));
% else
%   %  n_samples=12;
%    % parpool(12);
% end
% 
% addpath(genpath('/home/zhw22003/YALMIP-master'))
% system('export PATH=$PATH:/home/zhw22003/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zhw22003/mosek'))
% addpath(genpath('/jet/home/zwei4/YALMIP-master'))
% system('export PATH=$PATH:/jet/home/zwei4/mosek/11.0/tools/platform/linux64x86/bin')
% addpath(genpath('/home/zwei4/mosek'))
syms y t t2 n;
ttt1=0; %TIME START POINT
N = 6;
a=1.2;
b=0;
t_optimal_b=[];
t_optimal=0;
%bb2=linspace(0,2,4);%
T=200;
dt=1;
t0=[0, 20, 40, 60, 80, 100];
t0=20;

[x, D, D2, D4,w] =finitediff(N,2);

% [D,~]=cheb(N);
%for bb1=1:length(bb2)
%b=bb2(bb1);
%w
% [yi,w] = clencurt(N);
h=eye(size(D));
h1=eye(size(D));
k2=(a.^2+b.^2)*h;
k=(a.^2+b.^2).^0.5*h1;
w1=w.^0.5;
W=diag(w1);
d=D+k;
WW=W*d;
% WW=WW(2:N,2:N);
% W=W(2:N,2:N);

%v
V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126

% D2=D*D;
% D2=D2(2:N,2:N);
% S = diag([0; 1 ./(1-yi(2:N).^2); 0]);
% D4 = (diag(1-yi.^2)*D^4 - 8*diag(yi)*D^3 - 12*D^2)*S;
% D4 = D4(2:N,2:N);
c=D2-k2;


%U
Cn=0;
NN=100;
k0 = 0.1;
Re= 500;
g=1-exp(-k0.*t); %Acc
g=exp(-k0.*t); %Dec
an = (pi*n).^2./Re;  %lately
m=diff(g,t);
a1=subs(m,t,0);
m2=diff(m,t);
a2=subs(m2,t,t2);
e=exp(an.*t2);
a3=simplify(e*a2);
a4=simplify(int(a3, t2, 0, t),'Steps',200);
a5=a4+a1;
aa=-2.*Re.*(-1).^n./(pi.*n).^3.*(a5)+Cn;
a6=simplify(aa.*exp(-an.*t).*sin(n.*y.*pi),'Steps',100);
b1=Re./6.*m.*(y.^3-y);
b2=g.*y;
bb=b1+b2;
a7 = 0;
for ii = 1:NN
    a7 = a7 + subs(a6, n, ii);
end
U=a7+bb;

mm=diff(U);
mm2=diff(mm);
% U_yi = (subs(U, y, yi(2:N)));
% m_yi = (subs(mm, y, yi(2:N)));
% m2_yi = (subs(mm2, y, yi(2:N)));

U_yi = (subs(U, y, x));
m_yi = (subs(mm, y, x));
m2_yi = (subs(mm2, y, x));
U_yi_function = matlabFunction(U_yi);
m_yi_function = matlabFunction(m_yi);
m2_yi_function = matlabFunction(m2_yi);


%speeding code
o=1/(1i*Re)*(D4-2*D2*k2+k2*k2);
oo=-inv(c);
ooo=(1/(1i*Re))*c;


for j3 = 1:length(t0)

    %A
    t_steps = (T-t0(j3))/dt;
    ttt = linspace(ttt1, ttt1+T-t0(j3), t_steps);
    A = zeros(2*(N-2), 2*(N-2), t_steps);%
    %U_value
    G = zeros(1, t_steps);
    t_steps_ttt1=ttt1/dt;
    AA = eye(size(V));
    for i = 1:t_steps  %step
        t_value=t0(j3)+i*dt;

        %debug to set t_value=0; 
        % t_value=0;

%         U_vec = U_yi_function(t_value);
%         Uy_vec = m_yi_function(t_value);
%         Uyy_vec = m2_yi_function(t_value);

        U_yi_value = diag(U_yi_function(t_value));% t from t0 to t0+T
        m_yi_value = diag(m_yi_function(t_value));
        m2_yi_value = diag(m2_yi_function(t_value));

        %L
        Los=oo*(o-a*U_yi_value*c+a*m2_yi_value);%!
        Lsq=a*U_yi_value-ooo;
        Lc=b*m_yi_value;%diag
        L_value=[Los,zeros(N-2,N-2);Lc,Lsq];
        L=V*(-1i*L_value)/V;%new added



        A(:, :, i) = L; 
        AA = expm(-1i*L_value*dt)*AA;
        gh1=double(V*AA/V);
        G(i)=norm(gh1)^2;
    end

    G_storage{j3} = G;
    max_trans=max(G);
    %-----given tt value 
    yalmip('clear');
    I=eye(2*(N-2));
    P = cell(t_steps, 1);
    constraints=[];
    tt=sdpvar(1);

    for j1=1:t_steps
        P{j1}=sdpvar(2*(N-2),2*(N-2),'hermitian','complex');
        constraints = [constraints,I<=P{j1}<=tt*I];
    end

    scaling=Re;
    E=I*scaling;
    A_new=A*scaling;

    for j1=1:t_steps-1
        dP_dt = (P{j1 + 1}-P{j1})/dt;
        % constraints = [constraints, A(:,:,j1)'*P{j1}*E(:,:,j1)+E(:,:,j1)'*P{j1}* A_new(:,:,j1) + E(:,:,j1)'*dP_dt*E(:,:,j1) <=0];
        constraints = [constraints, A_new(:,:,j1)'*P{j1}*E+E'*P{j1}* A_new(:,:,j1) + E'*dP_dt*E <=0];

    end

    sdp_option=sdpsettings('solver','mosek');
    constraints=[constraints];
    % sdp_option.mosek.MSK_DPAR_ANA_SOL_INFEAS_TOL=1e-10;
    diagnostics = optimize(constraints,tt,sdp_option);


    for j1=1:t_steps
        eig_min(j1)=min(eig(value(P{1})));
        eig_max(j1)=max(eig(value(P{j1})));
    end


    if diagnostics.problem == 0
        % SUCCESS
        P_optimal = value(P);%eig (P)
        t_optimal = value(tt);
        disp('All optimal tt values:');
        disp(t_optimal);
        all_t_optimal(j3) = t_optimal;
        all_P_optimal{j3} = P_optimal;

    % Combine A matrices (t0 to T)
        A_full = cat(3, L, A);  % [A0, A1, A2, ..., A_180]
        num_trajectories = 3;

% Preallocate storage for V_lyap results
        all_V_lyap = zeros(num_trajectories, t_steps+1);

     for traj_idx = 1:num_trajectories
    % Random initial condition (complex)
        n = size(L,1);
        x0 = randn(n,1) + 1i*randn(n,1);
        x0 = x0 / norm(x0);  % Normalize
    
        % Simulate state evolution
        x = zeros(n, t_steps+1);
        x(:,1) = x0;
        for k = 1:t_steps
            x(:,k+1) = expm(A_full(:,:,k)*dt) * x(:,k);
        end
    
        % Time vectors
        t_sim = (0:dt:(t_steps*dt));  % Simulation times [t0, t0+T]
        t_P = linspace(0, T-t0, t_steps);  % P-matrix times

        % Compute x'Px at each time
        V_lyap = zeros(1, t_steps+1);
        for k = 1:(t_steps+1)
            if k == 1
                % Use first P matrix at t0
                P_mat = all_P_optimal{j3}{1};
            else
                % Find nearest P matrix in time
                [~, idx] = min(abs(t_P - t_sim(k)));
                P_mat = all_P_optimal{j3}{idx};
            end
            V_lyap(k) = real(x(:,k)' * P_mat * x(:,k));
        end
        all_V_lyap(traj_idx, :) = V_lyap;
     end
  
        % Compute eigenvalues/vectors for this t0
        eigvals = cell(t_steps, 1);
        eigvecs = cell(t_steps, 1);
        for j1 = 1:t_steps
            P_numeric = double(P_optimal{j1});
            [V_eig, D_eig] = eig(P_numeric);
            eigvals{j1} = diag(D_eig);
            eigvecs{j1} = V_eig;

        end
        
        all_eigenvalues{j3} = eigvals;
        all_eigenvectors{j3} = eigvecs;

       if ~isempty(all_eigenvalues{j3})
        %         time_points = t0(j3) + (1:length(all_eigenvalues{j3}))*dt;
                time_points = (1:length(all_eigenvalues{j3}))*dt;
                min_eig = cellfun(@min, all_eigenvalues{j3});
                max_eig = cellfun(@max, all_eigenvalues{j3});
        end
      
        bound = max_eig'.*G;
%         figure
%         hold on;
%         semilogy(t_sim, V_lyap, '--', 'LineWidth', 1.5, 'Color', 'r');
%         semilogy(time_points, bound, '-', 'LineWidth', 1.5, 'Color', 'b');
%         xlim([0 200]);
%         set(gca, 'YScale', 'log');
%         legend('show', 'Location', 'best');
        
    else
        % FAILED
        all_t_optimal(j3) = NaN;
        all_P_optimal{j3} = {};
        all_eigenvalues{j3} = {};
        all_eigenvectors{j3} = {};
        disp('Optimization failed.');
        disp(diagnostics.info);
    end

    t_plot=t_optimal*ones(size(ttt));
%     hold on
%     plot(ttt,t_plot,'--');
%     %axis control
%     xlim([0 200]);
%     set(gca, 'YScale', 'log');
%     ylim([10^-1 10^7]);% ylim([10^-1 10^5]) for fig7
%     set(gca, 'YTick', [10^-1 10^1 10^3 10^5 10^7]);%set(gca, 'YTick', [10^-1 10^1 10^3 10^5]); for fig7
end


% -------------------------------------------------------------------------

% figure;
% for j3 = 1:length(t0)
%     % Calculate the number of time steps for this segment
%     t_steps = (T - t0(j3)) / dt;
%     
%     % Generate time vector from 0 to T - t0(j3)
%     x_vals = linspace(0, T - t0(j3), t_steps);
%     
%     % Create constant line at the optimal tt value for this segment
%     y_vals = all_t_optimal(j3) * ones(1, t_steps);
%     
%     % Plot the horizontal segment
%     plot(x_vals, y_vals, '--');
%     hold on;
% end
% 
% xlim([0 200]);
% set(gca, 'YScale', 'log');
% ylim([10^-1 10^7]);
% set(gca, 'YTick', [10^-1 10^1 10^3 10^5 10^7]);
% xlabel('Time');
% ylabel('Optimal tt');
% title('Optimal tt over Time Windows');
% grid on;
% hold off;
% 
save('data_N_32.mat');
% save('P_matrices.mat', 'all_P_optimal', 't0');
% saveas(gcf, 'myplot.png');

%% Eigenvalue over time
% figure;
% hold on;
% cmap = parula(length(t0));  % Colormap with 9 colors (one per t0 value)
% 
% % Plot each curve with color based on t0 value
% for j3 = 1:length(t0)
%     if ~isempty(all_eigenvalues{j3})
%         time_points = t0(j3) + (1:length(all_eigenvalues{j3}))*dt;
%         min_eig = cellfun(@min, all_eigenvalues{j3});
%         max_eig = cellfun(@max, all_eigenvalues{j3});
%         semilogy(time_points, min_eig, 'Color', cmap(j3,:), 'LineWidth', 1.5);
%         semilogy(time_points, max_eig, '--', 'Color', cmap(j3,:), 'LineWidth', 1.5);
%     end
% end
% 
% % Configure colorbar with exact t0 values
% c = colorbar;
% colormap(cmap);
% 
% % Set discrete ticks at centers of each color segment
% c.Ticks = linspace(0, 1, length(t0));  % 9 ticks from 0 to 1
% c.TickLabels = arrayfun(@(x) num2str(x), t0, 'UniformOutput', false);  % Exact t0 values
% 
% legend({'min \lambda(P)', 'max \lambda(P)'}, 'Location', 'best');
% 
% % Axis labels and formatting
% xlabel('t');
% ylabel('Eigenvalue of \bf{P}');
% grid on;
% box on;
% hold off;


% % Initialize storage
% all_eigenvalues = cell(length(t0), 1);
% all_eigenvectors = cell(length(t0), 1);
% 
% for j3 = 1:length(t0)
%     if isempty(all_P_optimal{j3})
%         all_eigenvalues{j3} = {};
%         all_eigenvectors{j3} = {};
%         continue;
%     end
%     
%     P_cell = all_P_optimal{j3};
%     t_steps = length(P_cell);
%     eigvals = cell(t_steps, 1);
%     eigvecs = cell(t_steps, 1);
%     
%     for j1 = 1:t_steps
%         P_numeric = double(P_cell{j1});
%         [V_eig, D_eig] = eig(P_numeric);
%         eigvals{j1} = diag(D_eig);
%         eigvecs{j1} = V_eig;
%     end
%     
%     all_eigenvalues{j3} = eigvals;
%     all_eigenvectors{j3} = eigvecs;
% end
% 
% % Save results
% save('eigen_results.mat', 'all_eigenvalues', 'all_eigenvectors');

%% eigenvector contour plot
% for j3 = 1:length(t0)
%     if ~isempty(all_eigenvectors{j3})
%         % Full spatial grid (including boundaries)
% %         x_full = [-1,x,1];  % Size N x 1
%         Nx = 32;
% %          x_full = linspace(0,2*pi/a,Nx);
%         kx = a; kz = b;
%         k2 = kx^2 + kz^2;
%         
%         % Choose the first time step for mode visualization
%         j1 = 1;
%         t_current = t0(j3) + j1*dt;
%         
%        % After getting eigenvector q
%         [~, min_idx] = min(all_eigenvalues{j3}{j1});
%         q = all_eigenvectors{j3}{j1}(:, min_idx);
%         
%         % Split state vector
%         v_interior = q(1:(N-2));
%         wy_interior = q((N-2)+1:end);
%         
%         % Pad with boundary zeros (now size N)
%         v_full = [0; v_interior; 0];
%         wy_full = [0; wy_interior; 0];
%         
%         % Compute ∂v/∂y at INTERIOR points
%         Dv_interior = D * v_interior;  % D: (N-2)x(N-2) matrix
%         Dv_full = [0; Dv_interior; 0]; % Pad boundaries
%         
%         % Compute velocity components (now all size N)
% %         kx = a; kz = b;
% %         k2 = kx^2 + kz^2;
%         u_hat = (1/k2) * (1i*kx * Dv_full - 1i*kz * wy_full);
%         v_hat = v_full;
%         w_hat = (1/k2) * (1i*kz * Dv_full + 1i*kx * wy_full);
%         
% %         u = real(u_hat.*exp(1i*kx*x_full));
% 
%         y_grid = linspace(-1, 1, N)';
%         x_grid = linspace(0, 2*pi/a, Nx)';
%         [X, Y] = meshgrid(x_grid, y_grid);
%         
%         % Reconstruct physical space perturbation
%         u_physical = real(u_hat .* exp(1i*a*X));
%         v_physical = real(v_hat .* exp(1i*a*X));
%         w_physical = imag(w_hat .* exp(1i*a*X));
%         max_abs = max(abs(u_physical(:)));
%         c_limits = [-max_abs, max_abs];
%         % Plot
%         figure(j3);
%         contourf(X, Y, u_physical, 20, 'LineColor', 'none');
%         colorbar;
%         colormap(bluewhitered);
%         caxis(c_limits);
%         xlabel('Streamwise position (x)');
%         ylabel('Wall-normal position (y)');
%         title(sprintf('Perturbation u at t=%.1f (\\alpha=%.2f, \\beta=%.2f)', t0, a, b));
%         set(gca, 'FontSize', 12);
%         
%        
%     end
% end
% 
% 

% function [x, D1, D2, D4,w] =finitediff(N,L)
% 
% %Differentiation matrix using finite difference scheme.
% %This is suitable for Dirichlet boundary condition v(x=L/2)=v(x=-L/2)=0 at
% %the boundary and Neumann boundary condition v'(x=L/2)=v'(x=-L/2)=0. 
% 
% dx=L/N;%get grid spacing. 
% 
% x=linspace(-L/2,L/2,N);%get the grid point location
% 
% x=x(2:end-1); %delete the first and last points that are zero. 
% 
% w=dx*ones(1,N-2);%integration weighting 
% 
% N_diff=N-2;%The size of differentiation matrices
% 
% %First order derivative based on central difference
% %f'(x_i)=(f(x_{i+1})-f(x_{i-1})/(2*dx)
% %We also use the boundary condition that x_0=x_N=0
% D1_stencil=diag(-1*ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1);
% D1=D1_stencil/(2*dx);
% 
% %Second order derivative based on central difference
% %f''(x_i)=(f(x_{i+1})-2f(x_i)+f(x_{i-1})/(dx^2)
% %We also use the boundary condition that x_0=x_N=0
% D2_stencil=(diag(-2*ones(1,N_diff))+diag(ones(1,N_diff-1),1)+diag(ones(1,N_diff-1),-1));
% D2=D2_stencil/dx^2;
% 
% %Forth order derivative based on central difference
% %f''''(x_i)=(f(x_{i+2})-4f(x_{i+1})+6f(x_i)-4f(x_{i-1})+f(x_{i-2})/(dx^4)
% %This differentiation matrix only go through x_1 up to x_{N-1}
% D4_stencil=(diag(6*ones(1,N_diff))+diag(-4*ones(1,N_diff-1),1)+diag(-4*ones(1,N_diff-1),-1)...
%     + diag(ones(1,N_diff-2),2)+diag(ones(1,N_diff-2),-2));
% 
% %Here, we use the Neumann boundary condition that x_0'=x_N'=0
% %such that x_1=x_{-1} and x_{N-1}=x_{N+1} for the ghost points. Then we
% %also use the condition x_0 and x_N=0 to express all values based on x_1 up
% %to x_{N-1}
% D4_stencil(1,1)=D4_stencil(1,1)+1;
% D4_stencil(end,end)=D4_stencil(end,end)+1;
% 
% D4=D4_stencil/dx^4;
% 
% end


% function newmap = bluewhitered(m)
% %BLUEWHITERED   Blue, white, and red color map.
% %   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
% %   to red colormap, with white corresponding to the CAXIS value closest
% %   to zero.  This colormap is most useful for images and surface plots
% %   with positive and negative values.  BLUEWHITERED, by itself, is the
% %   same length as the current colormap.
% %
% %   Examples:
% %   ------------------------------
% %   figure
% %   imagesc(peaks(250));
% %   colormap(bluewhitered(256)), colorbar
% %
% %   figure
% %   imagesc(peaks(250), [0 8])
% %   colormap(bluewhitered), colorbar
% %
% %   figure
% %   imagesc(peaks(250), [-6 0])
% %   colormap(bluewhitered), colorbar
% %
% %   figure
% %   surf(peaks)
% %   colormap(bluewhitered)
% %   axis tight
% %
% %   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
% %   COLORMAP, RGBPLOT.
% 
% 
% if nargin < 1
%    m = size(get(gcf,'colormap'),1);
% end
% 
% 
% bottom = [0 0 0.5];
% botmiddle = [0 0.5 1];
% middle = [1 1 1];
% topmiddle = [1 0 0];
% top = [0.5 0 0];
% 
% % Find middle
% lims = get(gca, 'CLim');
% 
% % Find ratio of negative to positive
% if (lims(1) < 0) & (lims(2) > 0)
%     % It has both negative and positive
%     % Find ratio of negative to positive
%     ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
%     neglen = round(m*ratio);
%     poslen = m - neglen;
%     
%     % Just negative
%     new = [bottom; botmiddle; middle];
%     len = length(new);
%     oldsteps = linspace(0, 1, len);
%     newsteps = linspace(0, 1, neglen);
%     newmap1 = zeros(neglen, 3);
%     
%     for i=1:3
%         % Interpolate over RGB spaces of colormap
%         newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
%     end
%     
%     % Just positive
%     new = [middle; topmiddle; top];
%     len = length(new);
%     oldsteps = linspace(0, 1, len);
%     newsteps = linspace(0, 1, poslen);
%     newmap = zeros(poslen, 3);
%     
%     for i=1:3
%         % Interpolate over RGB spaces of colormap
%         newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
%     end
%     
%     % And put 'em together
%     newmap = [newmap1; newmap];
%     
% elseif lims(1) >= 0
%     % Just positive
%     new = [middle; topmiddle; top];
%     len = length(new);
%     oldsteps = linspace(0, 1, len);
%     newsteps = linspace(0, 1, m);
%     newmap = zeros(m, 3);
%     
%     for i=1:3
%         % Interpolate over RGB spaces of colormap
%         newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
%     end
%     
% else
%     % Just negative
%     new = [bottom; botmiddle; middle];
%     len = length(new);
%     oldsteps = linspace(0, 1, len);
%     newsteps = linspace(0, 1, m);
%     newmap = zeros(m, 3);
%     
%     for i=1:3
%         % Interpolate over RGB spaces of colormap
%         newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
%     end
%     
% end
% end
% 
% 
% function c = redblue(m)
% %REDBLUE    Shades of red and blue color map
% %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
% %   The colors begin with bright blue, range through shades of
% %   blue to white, and then through shades of red to bright red.
% %   REDBLUE, by itself, is the same length as the current figure's
% %   colormap. If no figure exists, MATLAB creates one.
% %
% %   For example, to reset the colormap of the current figure:
% %
% %             colormap(redblue)
% %
% %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
% %   COLORMAP, RGBPLOT.
% 
% %   Adam Auton, 9th October 2009
% 
% if nargin < 1, m = size(get(gcf,'colormap'),1); end
% 
% if (mod(m,2) == 0)
%     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
%     m1 = m*0.5;
%     r = (0:m1-1)'/max(m1-1,1);
%     g = r;
%     r = [r; ones(m1,1)];
%     g = [g; flipud(g)];
%     b = flipud(r);
% else
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
% 
% c = [r g b]; 
% 
% end
