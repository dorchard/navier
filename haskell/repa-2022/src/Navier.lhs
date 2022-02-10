> {-# LANGUAGE NoMonomorphismRestriction #-}

> module Navier where

> import Navier.Boundary
> import Navier.Datadefs
> import Navier.Init
> import Navier.Params
> import Navier.Output
> import Navier.Simulation

> import Data.Array.Repa

> import System.Environment
> import Debug.Trace
> import Text.Printf

> main = do
>   argv <- getArgs
>   let output_name = if ((length argv) >= 1) then argv!!0 else "foo"
>       output_freq = if ((length argv) >= 2) then read $ argv!!1 else 1
>       dims  = Z :. (jmax+2) :. (imax+2)

>       -- initialise arrays
>       u = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) ui)
>       v = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) vi)
>       f = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) 0.0)
>       g = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) 0.0)
>       rhs = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) 0.0)
>       p = fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) 0.0)

>       -- initialise the "flag" array, marking a circular peg 
>       (flag, ibound) = init_flag imax jmax delx dely
>       ifluid = (imax * jmax) - ibound

>       -- initial boundary conditions on the velocity components
>       (u', v') = apply_boundary_conditions u v flag imax jmax ui vi      
>       output = mainLoop 0.0 0 del_t 0.0 u' v' p f g rhs flag ifluid ibound output_freq

>   -- write frames out
>   write_ppm u v p flag imax jmax delx dely output_name 0 output_freq
>   mapM (\(u, v, p, i) -> 
>          write_ppm u v p flag imax jmax delx dely output_name (i+1) output_freq) output

> mainLoop :: Double -> Int -> Double -> Double -> Array DIM2 Double -> Array DIM2 Double 
>             -> Array DIM2 Double -> Array DIM2 Double -> Array DIM2 Double -> Array DIM2 Double
>             -> Array DIM2 Int -> Int -> Int -> Int
>             -> [(Array DIM2 Double, Array DIM2 Double, Array DIM2 Double, Int)]
> mainLoop t iters del_t res u v p f g rhs flag ifluid ibound output_freq
>            | t < t_end = 
>              let
>                del_t' = set_timestamp_interval del_t u v
>                (f', g') = compute_tentative_velocity u v f g flag del_t'
>                rhs' = compute_rhs f' g' rhs flag del_t'

>                -- the heavy computational work mostly happens in poisson
>                (p', res', itersor) = poisson p rhs' flag ifluid res

>                (u', v') = update_velocity u v f' g' p' flag del_t'
>                (u'', v'') = apply_boundary_conditions u' v' flag imax jmax ui vi

>                -- output message
>                msg = printf "%d t:%.8f del_t:%.8f, SOR iters:%3d, res:%.8f, bcells:%d"
>                        iters (t+del_t') del_t' itersor res' ibound           
>                frameOutput = (u'', v'', p', iters)
>              in
>                trace msg $
>                if (iters `mod` output_freq == 0) then
>                    frameOutput:(mainLoop (t+del_t') (iters+1) del_t' res' u'' v'' p' f' g' rhs' flag ifluid ibound output_freq)
>                else
>                    mainLoop (t+del_t') (iters+1) del_t' res' u'' v'' p' f' g' rhs' flag ifluid ibound output_freq
>            | otherwise = [(u, v, p, iters)]
