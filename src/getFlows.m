function [varargout] = getFlows(varargin)
% getFlows is a Sim3Tanks function. This function returns a data table with
% the values of the flow variables.

% Written by Arllem Farias, February/2024.
% Last update February/2024 by Arllem Farias.

%==========================================================================

if(nargin()<1)
    error(errorMessage(01));
elseif(nargin()>1)
    error(errorMessage(02));
end

if(isa(varargin{1},'Sim3TanksClass'))
    objSim3Tanks = varargin{1};
else
    error(errorMessage(07));
end

%==========================================================================
global SIM3TANKS_LISTS; %#ok<*GVMIS>

if(isempty(SIM3TANKS_LISTS))
    error(errorMessage(04));
else
    LIST_OF_FLOWS = SIM3TANKS_LISTS.LIST_OF_FLOWS;
end
%==========================================================================

q = objSim3Tanks.getFlowVariables();
if(~isempty(q))
    q = array2table(q,'VariableNames',LIST_OF_FLOWS);
end

varargout{1} = q;

end