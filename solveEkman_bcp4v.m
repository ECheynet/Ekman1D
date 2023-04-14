function [sol4c] = solveEkman_bcp4v(f,K,ug,z,opts,varargin)

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('Kp',zeros(size(K)));
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
Kp = p.Results.Kp ; % derivative of eddy diffusibity for momentum


solinit = bvpinit(z, [ug; 0;0;0]);
sol4c = bvp4c(@bvpfcn, @bcfcn, solinit, opts);



%% Nested functions
    function res = bcfcn(ya,yb)
        %         res = [ya(1);yb(1)-ug;...
        %             ya(2);yb(2);...
        %             ya(3);yb(3);...
        %             ya(4);yb(4)];
        res = [ya(1);yb(1)-ug;...
            ya(2);yb(2)];
    end

    function dydz = bvpfcn(z1,y)
        
        % y1 = u
        % y2 = v
        % y3 = u'
        % y4 = v'
        
        K1= interp1(z,K,z1);
        K1_p= interp1(z,Kp,z1);
        
        dydz = zeros(4,1);
        dydz(1) = y(3);
        dydz(2)= y(4);        
        dydz(3) = -1/K1.*(f.*y(2) + K1_p.*y(3));
        dydz(4) = -1/K1.*(f.*(ug-y(1)) + K1_p.*y(4) );
        
        
        
        
    end


end

