> module Navier.Params where

Boundary flags

> xlength = 22.0       -- Width of simulated domain
> ylength = 4.1        -- Height of simulated domain
> imax = 660 ::Int     -- Number of cells horizontally
> jmax = 120 :: Int    -- Number of cells vertically

> delx = xlength/660.0
> dely = ylength/120.0

> t_end = 40.0         -- Simultion runtime
> del_t = 0.003        -- (default) duration of each timestep
> tau = 0.5            -- safety factor for timestep control

> itermax = 10::Int   -- max number of SOR iterations
> eps = 0.001          -- stopping error threshold for SOR
> omega = 1.7          -- relaxation parameter for SOR
> gamma = 0.9          -- upwind differencing factor in PDE discretisation

> re = 150.0           -- Reynolds number
> ui = 1.0             -- Initial X velocity (West-east flow)
> vi = 0.0             -- Initial Y velocity

> intCast :: Int -> Double
> intCast = fromInteger . toInteger

>                               