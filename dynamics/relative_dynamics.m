function [Ad, Bd, state] = relative_dynamics(dt, config)
%RELATIVE_DYNAMICS  Unified interface for relative-motion prediction models.
%
%   [Ad, Bd, state] = relative_dynamics(dt, config)
%
%   Dispatches to either CWH (circular) or YA (elliptic) prediction model
%   based on config.dynamics_mode.
%
%   Inputs:
%       dt     - time step [s]
%       config - struct with fields:
%           .dynamics_mode  - 'CWH' or 'YA'
%
%         For CWH mode:
%           .n_orbit        - mean motion [rad/s]
%
%         For YA mode:
%           .orbit_cfg      - struct from elliptic_orbit_config()
%           .nu             - current true anomaly [rad]
%
%   Outputs:
%       Ad    - 6x6 state transition matrix
%       Bd    - 6x3 input matrix
%       state - updated state struct (pass back to next call)
%           For CWH: state = config (unchanged)
%           For YA:  state = config with .nu updated to end-of-step value

switch upper(config.dynamics_mode)
    case 'CWH'
        [Ad, Bd] = cwh_stm(dt, config.n_orbit);
        state = config;

    case 'YA'
        [Ad, Bd, nu1] = ya_stm(dt, config.orbit_cfg, config.nu);
        state = config;
        state.nu = nu1;  % advance true anomaly for next call

    otherwise
        error('relative_dynamics:unknownMode', ...
              'Unknown dynamics_mode "%s". Use "CWH" or "YA".', ...
              config.dynamics_mode);
end
end
