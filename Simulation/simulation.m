clear
syms y t t2 n;
N = 32;
a=1.2;b=0; 
% a=0;b=1.6;
T=200;
dt=1;
t0=[0, 10, 20, 40, 60, 80, 100];

[x, D, D2, D4,w] =finitediff(N,2);

h=eye(size(D));
h1=eye(size(D));
k2=(a.^2+b.^2)*h;
k=(a.^2+b.^2).^0.5*h1;
w1=w.^0.5;
W=diag(w1);
d=D+k;
WW=W*d;

V=[WW,zeros(N-2,N-2);zeros(N-2,N-2),W]; %126*126
c=D2-k2;

% Laminar flow
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
max_eigA = cell(1, length(t0)); % Cell array to store max eigenvalues for each t0
min_eigA = cell(1, length(t0)); % Cell array to store max eigenvalues for each t0
%A
for  j = 1:length(t0)
    t_steps = (T-t0(j))/dt;
    ttt = linspace(0, T-t0(j), t_steps);
    G = zeros(1, t_steps);
    AA = eye(size(V));
    max_eig_current = zeros(1, t_steps);

    %U_value
    for i = 1:t_steps  %step
        t_value=t0(j)+i*dt;          
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
        eigA = eig(L);
        
        % Store maximum eigenvalue magnitude
        max_eig_current(i) = max(real(eigA));
%         L_value=double(L);
        AA =expm(-1i*L_value*dt)*AA;


        %G
     gh1=double(V*AA/V);
        G(i)=norm(gh1)^2; 
 
    end
    G_storage{j} = G;
    max_trans=max(G);
    max_eigA{j} = max_eig_current;
%     min_eigA{j} = min_eig_current;
    plot(ttt,G,'-');
    hold on
end

%axis control
xlim([0 200]);
set(gca, 'YScale', 'log');
ylim([10^-1 10^7]);
set(gca, 'YTick', [10^-1 10^1 10^3 10^5 10^7]);
