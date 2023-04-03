function [sol4c,Km2,Kh2] = scm_bcp4v(latitude,para,z,opts)

%% Get parameters
kappa = 0.4;
Omega = 7.29e-5;
f = 2*Omega*sind(abs(latitude));
h = para.h;
Km = para.Km;
Kh = para.Kh;
alpha = para.alpha;
model = para.model;
L = para.L; % Obukhov length
u_star = para.u_star; % Obukhov length
bc_u = para.bc_u(:);
bc_v = para.bc_v(:);
bc_theta = para.bc_theta(:);
ug = bc_u(1);
if ~isempty(bc_theta),
    theta_top = bc_theta(1);
    theta_bottom = bc_theta(2);
end
%% Parametrize Km and Kh if necessary
if isempty(Km)||isempty(Kh)
    [Km,Kh,Kh_p,Km_p] = similarityFun(z,L,kappa,u_star,alpha,model);
end
%% Initial conditions
solinit = bvpinit(z, [bc_u(:);bc_v(:);bc_theta(:)]);

if ~isempty(bc_theta)
    sol4c = bvp4c(@bvpfcn_w_theta, @bcfcn_w_theta, solinit, opts);
else
    sol4c = bvp4c(@bvpfcn_wo_theta, @bcfcn_wo_theta, solinit, opts);
end

Km2 = interp1(z,Km,sol4c.x);
Kh2 = interp1(z,Kh,sol4c.x);

%% Nested functions
    function res = bcfcn_w_theta(y_bottom,y_top)
        res = [y_bottom(1);y_top(1)-ug;...
            y_bottom(2);y_top(2);...
            y_bottom(3)-theta_bottom(1);y_top(3)-theta_top];
    end
    function res = bcfcn_wo_theta(y_bottom,y_top)
        res = [y_bottom(1);y_top(1)-ug;...
            y_bottom(2);y_top(2)];
    end
    function dydz = bvpfcn_wo_theta(z1,y)
        % y1 = u        % y2 = v
        % y3 = u'        % y4 = v'
        Km1= interp1(z,Km,z1);
        dydz = zeros(4,1);
        dydz(1) = y(3);
        dydz(2)= y(4);
        dydz(3) = -1/Km1.*(f.*y(2));
        dydz(4) = -1/Km1.*(f.*(ug-y(1)));
    end
    function dydz = bvpfcn_w_theta(z1,y)
        
        % y1 = u        % y2 = v        % y3 = theta
        % y4 = u'        % y5 = v'        % y5 = theta'
        
        Km1= interp1(z,Km,z1);
        Kh1= interp1(z,Kh,z1);
        Km1_p= interp1(z,Km_p,z1);
        Kh1_p= interp1(z,Kh_p,z1);
        
        dydz = zeros(6,1);
        dydz(1) = y(4);
        dydz(2)= y(5);
        dydz(3)= y(6);
        dydz(4) = -1/Km1.*(f.*y(2) + Km1_p.*y(4));
        dydz(5) = -1/Km1.*(f.*(ug-y(1)) + Km1_p.*y(5));
        dydz(6)= -Kh1_p./Kh1.*y(6);
    end

    function [Km,Kh,Kh_p,Km_p] = similarityFun(z,L,kappa,u_star,alpha,model)
        
        Km_p = [];
        Kh_p = [];
        
        if isempty(model)
            if L<0
                phi_m =(1+15.2*abs(z/L)).^(-1/4); % Hogstrom (1988) from Dyer (1974)
                phi_h =0.95.*(1+11.6.*abs(z/L)).^(-1/2);% Hogstrom (1988)
                
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
                Kh_p =(20.*kappa.*u_star.*(1 - z/h).^(2.*alpha).*((58.*abs(z))/(5.*abs(L)) + 1).^(1/2))/19 - (40.*alpha.*kappa.*u_star.*z.*(1 - z/h).^(2.*alpha - 1).*((58.*abs(z))/(5.*abs(L)) + 1).^(1/2))/(19.*h) + (116.*kappa.*u_star.*z.*sign(z).*(1 - z/h).^(2.*alpha))/(19.*abs(L).*((58.*abs(z))/(5.*abs(L)) + 1).^(1/2));
                Km_p = kappa.*u_star.*(1 - z/h).^(2.*alpha).*((76.*abs(z))/(5.*abs(L)) + 1).^(1/4) - (2.*alpha.*kappa.*u_star.*z.*(1 - z/h).^(2.*alpha - 1).*((76.*abs(z))/(5.*abs(L)) + 1).^(1/4))/h + (19.*kappa.*u_star.*z.*sign(z).*(1 - z/h).^(2.*alpha))/(5.*abs(L).*((76.*abs(z))/(5.*abs(L)) + 1).^(3/4));
            elseif L>0
                phi_m =(1+4.8.*(z/L));% Hogstrom (1988)
                phi_h =0.95+4.5.*(z/L);% Hogstrom (1988)
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
                Kh_p = (kappa*u_star*(1 - z/h).^(2*alpha))./((9*z)./(2*L) + 19./20) - (9*kappa*u_star*z.*(1 - z./h).^(2*alpha))./(2*L.*((9*z)./(2*L) + 19./20).^2) - (2*alpha*kappa*u_star*z.*(1 - z./h).^(2*alpha - 1))./(h.*((9*z)./(2*L) + 19./20));
                Km_p = (kappa*u_star.*(1 - z./h).^(2*alpha))./((24*z)./(5*L) + 1) - (24*kappa*u_star*z.*(1 - z./h).^(2*alpha))./(5*L.*((24*z)./(5*L) + 1).^2) - (2*alpha*kappa*u_star*z.*(1 - z./h).^(2*alpha - 1))./(h.*((24*z)./(5*L) + 1));
            else
                phi_m = ones(size(z));
                phi_h = 0.95*ones(size(z));
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
                Km_p = kappa*u_star.*(1 - z./h).^(2*alpha) - (2*alpha*kappa*u_star*z.*(1 - z./h).^(2*alpha - 1))./h;
                Kh_p = (20*kappa*u_star.*(1 - z./h).^(2*alpha))./19 - (40*alpha*kappa*u_star*z.*(1 - z./h).^(2*alpha - 1))./(19*h);
            end
            
        elseif strcmpi(model,'YSU'),
            zeta = z./L;
            
            if L<0
                a=1;            b=2/3;            c=5;            d=0.35;
                phi_m  = 1 + a*zeta + b*zeta .*(1 + c - d*zeta) .*exp(-d*zeta);
                phi_h = phi_m;
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
            elseif L>0
                a=1;            b=2/3;            c=5;            d=0.35;
                phi_m  = 1 + a*zeta + b*zeta .*(1 + c - d*zeta)  .* exp(-d*zeta);
                phi_h = phi_m;
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
            else
                a=1;            b=2/3;            c=5;            d=0.35;
                phi_m  = 1 + a*zeta + b*zeta .*(1 + c - d*zeta) .* exp(-d*zeta);
                phi_h = phi_m;
                Km = kappa.*z.*u_star./phi_m .*(1-z/h).^(2*alpha);
                Kh = kappa.*z.*u_star./phi_h .*(1-z/h).^(2*alpha);
            end
        else
            error('Unknown model')
        end
    end

end

