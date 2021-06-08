function [zsea] = getz0_charnock(Uref,zref,alpha,varargin)
% [zsea] = getz0_charnock(Uref,zref,a_c) estimates the roughness of
% the sea surface knowing the mean wind speed Uref at a height zref. Note
% that the present function only works for near-neutral atmospheric
% stability.
% 
% Inputs
%   - Uref: scalar [1 x 1] : reference mean wind speed above the surface in
%   m/s at a height zref
%   - Uref: scalar [1 x 1] : reference height above the surface at which Uref is obtained. A typical
%   value of zref is 10 m.
%   - alpha: scalar [1 x 1] : Charnock coefficient (no dimension)
% Outputs
%   - zsea: scalar [1 x 1]: Estimate of the sear roughness in meters.
%   Typical values are between 1e-4 m and 1e-2 m.
% 
% Author: E. Cheynet - UiB, Norway - Last modified: 08-06-2021

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('kappa',0.41); % von karman constant
p.addOptional('g',9.81); % Acceleration of gravity
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa = p.Results.kappa ;
g = p.Results.g ; 

%%
gac = g./alpha;

% options = optimoptions('fsolve','Display','iter','PlotFcn',@optimplotfirstorderopt);
options = optimoptions('fsolve','Display','off');

fun = @charnock;
x0 = [2e-4];
zsea = fsolve(fun,x0,options);


if imag(zsea)~= 0
    warning('zsea is a complex-valued number. A value of 2e-4 m is used instead');
    zsea = 2e-4;
else
    fprintf('Charnock relationship leads to zsea = %1.1e m \n',zsea)
end




    

    function F = charnock(z0)
        F = z0.*gac.*(log(zref./z0)).^2 -(Uref.*kappa).^2;  
    end
end
