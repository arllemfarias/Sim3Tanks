function dxdt = sysDynamics(~,x,Qin,K,h0,Sc,Beta,w)
% sysDynamics is a Sim3Tanks function. This function describes the dynamic
% model of the three-tank system.
%
% Example:
%   dxdt = ysDynamics(t,x,Qin,K,h0,Sc,Beta,pNoise)
%       t   : time vector
%       x   : state vector
%       Qin : inlet flow rates
%       K   : states of the valves
%       h0  : transmission pipe height
%       Sc  : cross-sectional area of the tanks
%       Beta: constant value [mu*S*sqrt(2*g)]
%       w   : process noise vector

% https://github.com/e-controls/Sim3Tanks

%==========================================================================

if(nargin()<8)
    error(getMessage('ERR001'));
elseif(nargin()>8)
    error(getMessage('ERR002'));
end

%==========================================================================

Q1in = Qin(1);
Q2in = Qin(2);
Q3in = Qin(3);

Qx = sysFlowRates(x,K,h0,Beta);

Qa   = Qx(1);
Qb   = Qx(2);
Q13  = Qx(3);
Q23  = Qx(4);
Q1   = Qx(5);
Q2   = Qx(6);
Q3   = Qx(7);

dx1_dt = (1/Sc)*(Q1in - Qa - Q13 - Q1) + w(1);
dx2_dt = (1/Sc)*(Q2in - Qb - Q23 - Q2) + w(2);
dx3_dt = (1/Sc)*(Q3in + Q13 + Q23 + Qa + Qb - Q3) + w(3);

dxdt = [dx1_dt ; dx2_dt ; dx3_dt];

end