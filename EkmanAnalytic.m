function [u,v] = EkmanAnalytic(Ug,K,latitude,z,varargin)
%  [u,v] = EkmanAnalytic(Ug,K,latitude,z,varargin) computes the analytical
%  solutions of Ekman's equation for the lower atmosphere assuming a constant
%  eddy viscosity
%
% Inputs
%
% Ug: float [1x1]: Geostrophic wind velocity (m/s) component
% K:  float [1x1]: Eddy viscosity(m^2/s)
% latitude:  float [1x1]: latitude of the profile selected (in degrees)
% z:  float [1xNz]: height above ground (m)
%
% Outputs
% u:  float [1xNz]:  East-West wind velocity component (m/s)
% v:  float [1xNz]:  East-West wind velocity component (m/s)
%
% Author: E. Cheynet (UiB) - last modified 21-02-2020
%
% See also solveEkmanB.m


%% InputParser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Omega',7.29e-5); % rate of angular rotation of the Earth (rad/s)
p.parse(varargin{:});
Omega = p.Results.Omega; % rate of angular rotation of the Earth (rad/s)

%% Preallocation

f = 2*Omega*sind(abs(latitude));
d = sqrt(2*K./f);
gammaZ = z./d;
u = Ug.*(1-exp(-gammaZ).*cos(gammaZ));

if latitude >0
    v = Ug.*(exp(-gammaZ).*sin(gammaZ));
else
    v = -Ug.*(exp(-gammaZ).*sin(gammaZ));
end
end

