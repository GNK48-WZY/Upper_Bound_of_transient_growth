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

h=eye(size(D));
h1=eye(size(D));
k2=(a.^2+b.^2)*h;
k=(a.^2+b.^2).^0.5*h1;
w1=w.^0.5;
W=diag(w1);
d=D+k;
WW=W*d;

%v
V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126
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
end
