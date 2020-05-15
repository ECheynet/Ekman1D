function [u,v,ut,vt,t] = solveEkman(z,Ug,K,latitude,varargin)
%
% function [u,v,ut,vt,t] = solveEkman(z,Ug,K,K,latitude,varargin)
% numerically  solve Ekman's equation for the Ekman layer in the atmosphere
% using an explicit finite-difference scheme.
%
% Inputs
% z:  float [1xNz]: height above ground (m)
% Ug: float [1x1]: Geostrophic wind velocity (m/s) component
% K:  float [1x1] or [1xNz]: Eddy viscosity for the u component (m^2/s)
% K:  float [1x1] or [1xNz]: Eddy viscosity for the v component (m^2/s)
% latitude:  float [1x1]: latitude of the profile selected (in degrees)
%
% Optional inputs:
% 'Nmax':  integer [1x1]: maximal number of time-step for convergence study
% 'dt': float [1x1]: time step
% 'DT':  integer [1x1]: Number of steps between each the transient values of u and v are not
% stored. If DT = 100, the velocity data will be saved every 100 steps
% before convergence is reached.
% 'Omega':  float [1x1]: rate of angular rotation of the Earth (rad/s)
% ''critConv':  float [1x1]: change in m/s of u and v between each time
% step below which the steady-state is assumed reached.
%
% Outputs
% u:  float [1xNz]:  East-West wind velocity component (m/s)
% v:  float [1xNz]:  North-South wind velocity component (m/s)
% ut: float [NxNz]:  East-West wind velocity component stored every DT steps (m/s)
% vt: float [NxNz]:  North-South wind velocity component stored every DT steps (m/s)
% t:  float [NxNz]:  time stored every DT steps (m/s)
%
% Author: E. Cheynet (UiB) - last modified 21-02-2020
%
%  see also EkmanSpiral.m%




%% InputParser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Nmax',2e4);
p.addOptional('dt',[]);
p.addOptional('DT_save',100);
p.addOptional('Omega',7.29e-5);
p.addOptional('critConv',1e-6);
p.addOptional('method','Euler');
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
Nmax = p.Results.Nmax ; % Maximal number of time steps (for convergence)
dt = p.Results.dt ; % time step for RK4
Omega = p.Results.Omega; % rate of angular rotation of the Earth (rad/s)
critConv= p.Results.critConv; % convergence criterion
DT_save= p.Results.DT_save; % Intermediate value are stored every DT step
method= p.Results.method;
%% Preallocation and Initial conditions
Nz = numel(z);
f = 2*Omega*sind(latitude);
dz = [0,diff(z)];
ut = zeros(Nmax,Nz);
vt = zeros(Nmax,Nz);
e = ones(Nz, 1);
A = spdiags( [e -2*e e], -1:1, Nz, Nz);

if isempty(dt)
    dt = 0.4*dz.^2/max(K);
end

if (dt>0.48*dz.^2/max(K)), warning('Time step may be too large.'); end

%% Initial conditions
u = Ug.*ones(1,Nz);
v = zeros(1,Nz);

%% Main loop
myFun = @(Y,A,F) A*Y + F;
count= 1;
for time_step=1:Nmax % time step
    %  update u,v a nd time at each step
    oldU = u;
    oldV = v;
    
    if strcmpi(method,'Euler')
        [u,v] = Euler(u,v,dt,Ug,K,f,A,dz);
    elseif strcmpi(method,'RK4')
        [u,v] = solve_with_RK4(myFun,u,v,Ug,K,f,A,dt,dz);
    end
    
    % Boundary conditions
    [u,v] = applyBC(u,v,Ug,Nz);
    
    % Store data very N step
    if mod(time_step,DT_save)==0
        ut(count,:)=u;
        vt(count,:)=v;
        count = count+1;
        if time_step>1 && time_step<Nmax
            % Check convergence
            maxDiffU = max(abs(u-oldU));
            maxDiffV = max(abs(v-oldV));
            if maxDiffU < critConv && maxDiffV< critConv
                ut = ut(1:count-1,:);
                vt = vt(1:count-1,:);
                t = [0:count-2].*mean(dt);
                fprintf('The algorithm has converged \n')
                return
            end
        elseif time_step==Nmax
            ut = ut(1:count-1,:);
            vt = vt(1:count-1,:);
            t = [0:count-2].*mean(dt);
            warning('The algorithm did not converge')
        else
            error('Unknown error')
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [u,v] = Euler(u,v,dt,Ug,K,f,A,dz)
        
        % At each time step, compute the second derivative of u and v
        d2u = (A*u')'./dz.^2.*K;
        d2v = (A*v')'./dz.^2.*K;
        
        
        if numel(K)==1
            u = u + dt.*(d2u + f.*v) ;
            v = v + dt.*(d2v + f.*(Ug-u));
        else
            dK = [0,diff(K)];
            dU = [0,diff(u)];
            dV = [0,diff(v)];
            u = u + dt.*(d2u + f.*v + dK.*dU./dz.^2);
            v = v + dt.*(d2v + f.*(Ug-u)+ dK.*dV./dz.^2);
        end
    end
    function [u,v] = applyBC(u,v,Ug,Nz)
        % apply the upper boundary condition
        u(Nz) = Ug;
        v(Nz) = 0;
        
        % apply the no-slip lower boundary condition
        u(1) = 0;
        v(1) = 0;
    end
    function [u,v] = solve_with_RK4(Fun,u,v,Ug,K,f,A,dt,dz)
        
        
%         dt = repmat(w*(dz.^2)./K,2,1); % time step
        % At each time step, compute the second derivative of u and v
        d2u = K.*(A*u')'./dz.^2;
        d2v = K.*(A*v')'./dz.^2;
        
        if numel(K)==1
            F = [d2u; d2v + f.*Ug];
        else
            dK = [0,diff(K)];
            dU = [0,diff(u)];
            dV = [0,diff(v)];
            F = [d2u + dK.*dU./dz.^2; d2v + dK.*dV./dz.^2 + f.*Ug];
        end
        
        Y = [u;v];
        B = [0 f; -f 0];
        
        
        
        % Runge-Kutta of order 4
        k_1 = Fun(Y,B,F);
        k_2 = Fun(Y+0.5*dt.*k_1,B,F);
        k_3 = Fun(Y+0.5*dt.*k_2,B,F);
        k_4 = Fun(Y+k_3.*dt,B,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4).*dt;
        
        u = Y(1,:);
        v = Y(2,:);
%         dt = dt(1,:);
        
        
    end
end

