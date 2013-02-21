> {-# LANGUAGE QuasiQuotes #-}

> module Navier.Init where

> import Ypnos
> import Ypnos.Core.Grid

> import Data.Bits

> flagBoundary = [boundary| Double (*i, -1) -> _cb -- Mark N & S boundary cells
>                                   (*i, +1) -> _cb
>                                   (-1, *j) -> _cb -- Mark E & W boundary cells
>                                   (+1, *j) -> _cb |]

> flagBoundaryCells = [fun| X*Y:| _  n _ |  
>                               | e @c w |
>                               | _  s _ | ->  if (c .&. _cf == 0) then 
>                                               Reduce 1 
>                                              (if (e .&. _cf) then c .|. _bw else
>                                               if (w .&. _cf) then c .|. _be else
>                                               if (n .&. _cf) then c .|. _bs else
>                                               if (s .&. _cf) then c .|. _bn else c)
>                                              else Reduce 0 c |]

> init_flag flag imax jmax delx dely =
>     let
>         mx = 20.0/41.0*jmax*dely
>         my = mx
>         init g = let (i, j) = cursor g
>                      x = (i-0.5)*delx - mx
>                      y = (j-0.5)*dely - my
>                  in if (x*x + y*y <= rad1*rad1)
>                     then _cb 
>                     else _cf 
>     in 
>        (runReduce flagBoundaryCells) . (addBoundary flagBoundary) . (run init) $ flag