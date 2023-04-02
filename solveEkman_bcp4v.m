function [sol4c] = solveEkman_bcp4v(f,K,ug,z,opts)




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
        
        dydz = zeros(4,1);
        dydz(1) = y(3);
        dydz(2)= y(4);
        dydz(3) = -1/K1.*f.*y(2);
        dydz(4) = -1/K1.*f.*(ug-y(1));
        
        
        
        
        
        
    end


end

