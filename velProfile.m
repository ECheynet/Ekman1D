
function [U] = velProfile(u_star,z,z0,L,varargin)
%  [U] = velProfile(u_star,z,z0,L,varargin) computes the wind speed
%  profile in the surface layer under neutral and non-neutral thermal
%  stratification, providing that z/L>-2 and z/L<1. 
% The profile is computed by direct numerical integration of the
% non-dimensional similarity profile in the surface layer.
% 
% INPUTS:
% u_star: scalar [1x1]: friction velocity (m/s)
% z: vector [1xNz] or [Nzx1]: Height above the surface (m)
% z0: scalar [1x1]: roughness length (m)
% L: scalar [1x1]: Obukhov length (m)
% varargin
% kappa: scalar [1x1]: von Karman constant, by default equal to 0.40
% C: vector [1x2] constant parameters used in the non-dimensional velocity profile
% 
% 
% OUTPUT
% U: vector [1xNz] or [Nzx1]: Mean wind speed profile
% 
% Author: E. Cheynet - UiB - last modified: 17.12.2019
% 
% see also tempProfile
%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('kappa',0.40);
% p.addOptional('C',[5 15.2]); % from Dyer
p.addOptional('C',[5 16]); % from Dyer
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = p.Results.kappa;
C = p.Results.C;
if abs(nanmean(z./L))<1e-3 % near neutral conditions
    U = zeros(1,numel(z));
    U(z>z0) = u_star./kappa.*log(z(z>z0)./z0);
    U(z<=z0) = u_star./kappa.*log(z0./z0);
    
else
    
    Z = logspace(log10(0.99*z0),log10(1.01*max(z)),2000);
    [phiM] = getPhiM(Z,L,C);
    U = u_star./kappa.*cumtrapz(Z,phiM./Z);
    U = interp1(Z,U,z);
end
    function [phi] = getPhiM(z,L,C)
        if L>=0
            phi = 1+C(1).*(z./L);
        else
            phi = (1+C(2).*abs(z./L)).^(-1/4);
        end
        
    end
end
