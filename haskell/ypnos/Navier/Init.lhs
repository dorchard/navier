> {-# LANGUAGE QuasiQuotes #-}

> module Navier.Init where

> import Ypnos
> import Ypnos.Core.Grid

> init_flag flag imax jmax delx dely =
>     let
>         init g = let (i, j) = cursor g
>                      x = (i-0.5)*delx - mx
>                      y = (j-0.5)*dely - my
>                  in if (x*x + y*y <= rad1*rad1) then _cb 
>                     else _cf 
>                         
>                         
>     in
>       
