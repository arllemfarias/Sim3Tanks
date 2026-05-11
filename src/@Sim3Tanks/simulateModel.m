function [ysol,xsol,qsol] = simulateModel(varargin)
% simulateModel is a Sim3Tanks method. This method simulates the dynamic
% behavior of the three-tank system defined by the user.
%
% The input parameter must follow the pair ('NAME',VALUE), where NAME must
% be Qp1, Qp2, Qp3, or Tspan, and VALUE must be a numeric type.
%
% If the pair ('Qp1',VALUE1) is omitted, then the declared value in the
% field objSim3Tanks.Model.PhysicalParam.PumpMaxFlow is used as default.
% The same is valid for the pairs ('Qp2',VALUE2) and ('Qp3',VALUE3). For
% the pair ('Tspan',VALUE), the default value is 0.1.
%
% NOTE: It is highly recommended to use simulation time increment as Tspan.
%
% The pair ('allSteps',true) enables the return of all intermediate steps
% of each simulation call. By default, its value is false, and only the
% final step is returned.
%
% The output arguments are the following vectors:
%   y = [h1,h2,h3,Q1in,Q2in,Q3in,Qa,Qb,Q13,Q23,Q1,Q2,Q3] : measurement vector
%   x = [h1,h2,h3] : state vector
%   q = [Q1in,Q2in,Q3in,Qa,Qb,Q13,Q23,Q1,Q2,Q3] : flow vector
%
% Examples of how to call the method:
%
%   tts = createSim3Tanks(); % create an object.
%   tts.setDefaultModel(); % configure the object to the default model.
%
%   [y,x,q] = tts.simulateModel(); % default values are used.
%
%   [y,x,q] = tts.simulateModel('Qp1',100); % only the value of Qp1 is
% % updated to 100.
%
%   [y,x,q] = tts.simulateModel('Qp2',110,'Tspan',0.2); % only the values
% % of Qp2 and Tspan are updated to 110 and 0.2, respectively.
%
%   [y,x,q] = tts.simulateModel('allSteps',true); % the variables y, x,
% % and q will have more than one line.
%
% See also createSim3Tanks, setDefaultModel, getSensorMeasurements,
%          getStateVariables, getFlowVariables, getFaultSignals,
%          getValveSignals.

% https://github.com/e-controls/Sim3Tanks

%==========================================================================

if(mod(nargin(),2) == 0)
    error(getMessage('ERR020'));
else
    objSim3Tanks = varargin{1};
end

%==========================================================================

[Param,ID] = checkPhysicalParam(objSim3Tanks);

Rtank = Param.(ID{1});
Hmax  = Param.(ID{2});
Rpipe = Param.(ID{3});
h0    = Param.(ID{4});
mu    = Param.(ID{5});
g     = Param.(ID{6});
Qmin  = Param.(ID{7});
Qmax  = Param.(ID{8});

Sc = pi()*(Rtank^2); % Cross-sectional area of the tanks (cm^2)
S  = pi()*(Rpipe^2); % Cross-sectional area of the pipes (cm^2)
Beta = mu*S*sqrt(2*g); % Constant value

%==========================================================================

% Default options
options.QP1 = Qmax;
options.QP2 = Qmax;
options.QP3 = Qmax;
options.TSPAN = 0.1;
options.ALLSTEPS = false;

% Check input options
for i = 2 : 2 : nargin()
    name = upper(varargin{i});
    if(isfield(options,name))
        options.(name) = varargin{i+1};
    else
        error(getMessage('ERR003'));
    end
end

% New options
Qp1 = options.QP1;
Qp2 = options.QP2;
Qp3 = options.QP3;
Tspan = options.TSPAN;
allSteps = options.ALLSTEPS;

%==========================================================================

opMode = checkOperationMode(objSim3Tanks);
[valveID,openingRate] = checkEnabledValves(objSim3Tanks);
[faultID,faultMag,offset] = checkEnabledFaults(objSim3Tanks);

K = zeros(size(opMode));

% ID = Sim3Tanks.LIST_OF_VALVES;
for i = 1 : numel(opMode)
    OP = opMode(i);
    EC = isempty(valveID{i});
    EF = isempty(faultID{i});

    % State machine to select the valve behavior
    if(~OP && EC && EF || OP && EC && EF)
        % fprintf('STATE 1 : %s = OP\n',ID{i});
        K(i) = opMode(i);

    elseif(~OP && EC && ~EF)
        % fprintf('STATE 2 : %s = f\n',ID{i});
        K(i) = faultMag(i);

    elseif(~OP && ~EC && EF || OP && ~EC && EF)
        % fprintf('STATE 3 : %s = K\n',ID{i});
        K(i) = openingRate(i);

    elseif(~OP && ~EC && ~EF || OP && ~EC && ~EF)
        % fprintf('STATE 4 : %s = K*(1-f)\n',ID{i});
        K(i) = openingRate(i)*(1-faultMag(i));

    elseif(OP && EC && ~EF)
        % fprintf('STATE 5 : %s = 1-f\n',ID{i});
        K(i) = 1-faultMag(i);

    else % Invalid state
        error(getMessage('ERR000'));
    end
end

%==========================================================================

Nx = numel(Sim3Tanks.LIST_OF_STATES);
Nq = numel(Sim3Tanks.LIST_OF_FLOWS);

Qin = [...
    K(1)*satSignal(Qp1,[Qmin Qmax]),...
    K(2)*satSignal(Qp2,[Qmin Qmax]),...
    K(3)*satSignal(Qp3,[Qmin Qmax]),...
    ];

x_vec = objSim3Tanks.getInternalStateVariables();

if(isempty(x_vec))

    x0 = objSim3Tanks.Model.InitialCondition;

    if(isempty(x0)||~isrow(x0)||numel(x0)~=numel(Sim3Tanks.LIST_OF_STATES))
        error(getMessage('ERR024'));
    end

    Qx = sysFlowRates(x0,K,h0,Beta);

    q = [Qin,Qx];

    [~,mNoise] = checkEnabledNoises(objSim3Tanks);

    y = sysMeasurements(x0,q,faultMag,offset,mNoise);

    objSim3Tanks.setInternalStateVariables(x0);
    objSim3Tanks.setInternalFlowVariables(q);
    objSim3Tanks.setInternalSensorMeasurements(y);
    objSim3Tanks.setInternalValveSignals(opMode');
    objSim3Tanks.setInternalFaultMagnitudes(faultMag');
    objSim3Tanks.setInternalFaultOffsets(offset(11:end)');
    objSim3Tanks.resetInternalSimulationTime();

else
    x0 = x_vec(end,:);
end

%==========================================================================

% Solver Configuration
options = odeset('MaxStep',Tspan,'RelTol',1e-6);

[pNoise,~] = checkEnabledNoises(objSim3Tanks);

model = @(t,x)sysDynamics(t,x,Qin,K,h0,Sc,Beta,pNoise);

[tsol,xsol] = ode45(model,[0 Tspan],x0,options);

if(all(isfinite(xsol),'all'))
    if(allSteps)
        xsol = satSignal(xsol(2:end,:),[0 Hmax]);
        tsol = tsol(2:end);
    else
        xsol = satSignal(xsol(end,:),[0 Hmax]);
        tsol = tsol(end);
    end
else
    error(getMessage('ERR007'));
end

numberOfSteps = size(xsol,1);
ysol = zeros(numberOfSteps,Nx+Nq);
qsol = zeros(numberOfSteps,Nq);

t0 = objSim3Tanks.getInternalSimulationTime(end);

for i = 1 : numberOfSteps

    Qx = sysFlowRates(xsol(i,:),K,h0,Beta);

    % Flows --> q = [Q1in,Q2in,Q3in,Qa,Qb,Q13,Q23,Q1,Q2,Q3]
    qsol(i,:) = [Qin,Qx];

    [~,mNoise] = checkEnabledNoises(objSim3Tanks);

    % Measurements --> y = [h1,h2,h3,Q1in,Q2in,Q3in,Qa,Qb,Q13,Q23,Q1,Q2,Q3]
    ysol(i,:) = sysMeasurements(xsol(i,:),qsol(i,:),faultMag,offset,mNoise);

    objSim3Tanks.pushInternalStateVariables(xsol(i,:));
    objSim3Tanks.pushInternalFlowVariables(qsol(i,:));
    objSim3Tanks.pushInternalSensorMeasurements(ysol(i,:));
    objSim3Tanks.pushInternalValveSignals(K');
    objSim3Tanks.pushInternalFaultMagnitudes(faultMag');
    objSim3Tanks.pushInternalFaultOffsets(offset(11:end)');
    objSim3Tanks.incrementInternalSimulationTime(t0+tsol(i));

end