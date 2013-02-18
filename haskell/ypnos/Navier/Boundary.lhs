> {-# LANGUAGE QuasiQuotes #-}
> {-# LANGUAGE TemplateHaskell #-}
> {-# LANGUAGE GADTs #-}

> {-# LANGUAGE NoMonomorphismRestriction #-}

> module Navier.Boundary where

> import Ypnos
> import Ypnos.Core.Grid

 import Navier

 apply_boundary_conditions :: Grid (Dim X :* Dim Y) b dyn Double ->
                              Grid (Dim X :* Dim Y) b dyn Dboule ->
                              Grid (Dim X :* Dim Y) b dyn Dboule ->

> xn = fst . size 
> yn = snd . size

> 
> ubound = [boundary| Double (-1, *j) g -> g!!!(0, j)           -- float from west
>                            (+1, *j) g -> g!!!((xn g) - 1, j)  -- float out to east
>                            (*i, -1) g -> g!!!(i, 0)
>                            (*i, +1) g -> g!!!(i, (yn g) - 1) |]

> vbound = [boundary| Double (-1, *j) g -> g!!!(0, j)
>                            (+1, *j) g -> g!!!((xn g) - 1, j)
>                            (*i, -1)   -> 0.0
>                            (*i, +1)   -> 0.0 |]

> no_slip_boundaries = [fun| X*Y:| _   t  _ | 
>                                | l  @c  r | 
>                                | _   b  _ | -> 0.0 |]

> apply_boundary_conditions u_data v_data flag imax jmax ui vi = 
>     let
>         u = listGrid (Dim X :* Dim Y) (0, 0) (imax, jmax) u_data ubound
>         v = listGrid (Dim X :* Dim Y) (0, 0) (imax, jmax) v_data vbound
>     in
>       (u, v)

