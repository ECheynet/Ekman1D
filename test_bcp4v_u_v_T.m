clearvars;close all;
h = 1000;
zmax = 0.99*h;
Nz = 200;
latitude = 53;
z0 = 0.001;
z = logspace(log10(z0),log10(zmax),Nz);

%% With explicit formulation for Km and temperature data
% Definitions of K1,K2, K3 and K4 following [1]
% anonymous function fo define K2,K3 and K4
getK_V1 = @(K0,K_hat,K_star,z_star,z)  K_hat + (K0-K_hat).*exp(z/z_star.*(log((K_star-K_hat)./(K0-K_hat))));
getK_V2 = @(a0,a1,K0,K_hat,K_star,z_star,z)  (a0 + a1.*z).*K0.*exp(z/z_star.*(log((K_star-K_hat)./(K0-K_hat))));

Nk = 4;
K = zeros(Nk,Nz);
K(1,:) = 1.3.*ones(1,Nz); % K1

K0 = 0.7; K_hat = 5.5; K_star = 4.5; z_star = 500;
K(2,:) = getK_V1(K0,K_hat,K_star,z_star,z); % K2

K0 = 1.4; K_hat = 0.7; K_star = 0.8; z_star = 500;
K(3,:) = getK_V1(K0,K_hat,K_star,z_star,z); % K3

a0 = 0.9; a1 = 2.7/50; K0 = 0.45; K_hat = 0.1; K_star = 0.3; z_star = 250;
K(4,:) = getK_V2(a0,a1,K0,K_hat,K_star,z_star,z); % K4

% z = linspace(1,h,Nz); % height vector
figure
for ii=1:Nk
    
    opts = bvpset('RelTol',0.01,'AbsTol',0.01,'Stats','on');
    
    clear para
    para.Km = K(ii,:) ; % No explicit formulation for Km --> MO theory is used
    para.Kh = para.Km ; % No explicit formulation for Km --> MO theory is used
    para.L = inf; % Obukhov length
    para.u_star = 0.25; % initial conditions
    para.h = h;
    para.alpha = 1;
    para.model = [];
    para.bc_theta = []; % boundary and initial conditions: [top-bottom]
    para.bc_u = [10 0];% boundary and initial conditions: [top-bottom]
    para.bc_v = [0 0]; % boundary and initial conditions: [top-bottom]
    
    [sol4c] = scm_bcp4v(latitude,para,z,opts);
    
    myZ = sol4c.x;
    myU = sol4c.y(1,:);

    plot(myU,myZ)
    hold on; box on;
    xlabel('u (m s^{-1})')
    ylabel('height (m)')
    set(gcf,'color','w')
end
legend('K1','K2','K3','K4','location','best')

%% Without explicit formulation for Km and temperature data
opts = bvpset('RelTol',0.01,'AbsTol',0.01,'Stats','on');

clear para
para.Km = [] ; % No explicit formulation for Km --> MO theory is used
para.Kh = para.Km ; % No explicit formulation for Km --> MO theory is used
para.model = [];
para.L = -1000; % Obukhov length
para.u_star = 0.25; % initial conditions
para.h = h;
para.alpha = 1;
para.bc_theta = [273 285]; % boundary and initial conditions: [top-bottom]
para.bc_u = [10 0];% boundary and initial conditions: [top-bottom]
para.bc_v = [0 0]; % boundary and initial conditions: [top-bottom]

%% Solve the single-column model
[sol4c,Km,Kh] = scm_bcp4v(latitude,para,z,opts);

% Process the outputs
if ~isempty(para.bc_theta),
    u = sol4c.y(1,:);
    v = sol4c.y(2,:);
    theta = sol4c.y(3,:);
    dudz = sol4c.y(4,:);
    dvdz = sol4c.y(5,:);
    dTdz = sol4c.y(6,:);
    wT = -Kh.*dTdz;
else
    u = sol4c.y(1,:);
    v = sol4c.y(2,:);
    dudz = sol4c.y(3,:);
    dvdz = sol4c.y(4,:);
end

uw = -Km.*dudz;
vw = -Km.*dvdz;
u_star = ((uw).^2 + (vw).^2).^0.25 ;

clf;close all
figure('position',[521   379   668   420]);
if ~isempty(para.bc_theta),
    tiledlayout(1,3,'TileSpacing','compact')
else
    tiledlayout(1,2,'TileSpacing','compact')
end

nexttile
meanU = sqrt(u.^2 +v.^2);
plot(meanU ,sol4c.x,'linewidth',1.2)
[~,indZ]=min(abs(sol4c.x-10));
hold on
uref_sim =  meanU(indZ);
u1 = velProfile(para.u_star,sol4c.x,z0,para.L);
uref_log =  u1(indZ);
plot(u1.*uref_sim./uref_log,sol4c.x,'k','linewidth',1.2);
ylabel('z (m)')
xlabel('$\overline{u}$ (m s$^{-1}$)','interpreter','latex')
grid on
legend('SCM','log profile')

nexttile
plot(uw,sol4c.x)
hold on
plot(vw,sol4c.x)
if ~isempty(para.bc_theta),
    plot(wT,sol4c.x,'linewidth',1.2)
    legend('uw','vw','wT','location','best')
else
    legend('uw','vw','location','best')
end
xlabel('(m$^2$ s$^{-2}$ or K m s$^{-1}$)','interpreter','latex')
grid on

if ~isempty(para.bc_theta),
    nexttile
    plot(theta,sol4c.x,'linewidth',1.2)
    xlabel('T (K)','interpreter','latex')
end
set(gcf,'color','w')
grid on

set(findall(gcf,'-property','FontSize'),'FontSize',12,'FontName','Times')