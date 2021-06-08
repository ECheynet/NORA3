function [newU,rmse] = interpU(oldZ,z0,oldU,newZ,u_star,varargin)
% [newU,rmse] = fitMeanU(z,z0,oldU,newZ,u_star) interpolates the mean wind speed
% profile measured or modelled in the atmospheric boundary layer (ABL) using the parabolic model of Deaves and Harris
% [1,2]. The interpolation is, therefore, non-linear.
% 
% Inputs
%   - oldZ: array [1 x N0] containig the height (in m) of original
%   data points.
%   - z0: scalar [1 x 1]: roughness length (in m) .
%   - oldU: array [1 x N0] containig the mean wind speed (in m/s) at the
%   heights oldZ
%   - newZ: array [1 x N1] containig the height (in m) of the interpolated
%   wind speed
%   - u_star: : scalar [1 x 1]: surface friction velocity (in m/s) .
% 
% Outputs
%   - newU: array [1 x N1] containing the mean wind speed (in m/s) at the
%   heights newZ
%   - rmse: scalar [1 x 1] containing the  root mean square error between fitted and measured mean wind
% speed profile.
% 
% References:
% 
% [1] Harris, R. I., & Deaves, D. M. (1980, November). The structure of strong winds, wind engineering in the eighties. In Proc. CIRIA Conf.
% [2] ESDU. (1985). ESDU 85020-Characteristics of atmospheric turbulence near the ground. Part II: single point data for strong winds (neutral atmosphere.
% 
% 
% Author: E. Cheynet - UiB, Norway - Last modified: 08-06-2021

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('kappa',0.41); % von karman constant
p.addOptional('h',1000); % height of the ABL
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = p.Results.kappa ;
h = p.Results.h ; 
% the h value does not affect the quality of the interpolation but the
% fitted coefficients. However, the function focuses on interpolating a
% wind speed profile not on identifying these coefficients in a reliable
% way.
%%
guess = getCoeff(h); % Get the first guess using the parameters from the Deaves and Harris model
myFun = @paraU;

% quick and dirty fit
if numel(oldZ)<numel(guess),
    warning('The problem may be over parametrized. You should use more data point for oldZ and oldU');
end
% 
warning off
C = nlinfit(oldZ,oldU,myFun,guess);
warning on

newU = myFun(C,newZ);

% Nearest neighbour extrapolation for newZ > max(oldZ);
newU(newZ>=oldZ(end))= oldU(end);

% Compute the root mean square error between fitted and measured mean wind
% speed profile.
rmse = RMSE(myFun(C,oldZ),oldU);

%%

    function A = paraU(coeff,z)
        % That is the Deaves and Harris wind speed profile
        A0 = u_star./kappa.*log(z./z0);
        z1 = z-z0;
        A1 = z1.*coeff(1);
        A2 = (z1.^2).*coeff(2);
        A3 = -(z1.^3).*coeff(3);
        A4 = -(z1.^4).*coeff(4);
        A = A0 + A1 +A2 +A3 + A4;
    end

    function C = getCoeff(h)
        % coefficients of the DEaves and harris model
        C(1) = 23/4./h;
        C(2) = -15/8/h.^2;
        C(3) = -4/3/(h).^3;
        C(4) = 1/4./(h).^4;
    end

    function [rmse] = RMSE(y1,y2)
        
        rmse = sqrt(nanmean((y1(:)-y2(:)).^2));
        
    end


end

