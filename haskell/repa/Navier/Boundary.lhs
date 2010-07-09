> {-# LANGUAGE NoMonomorphismRestriction #-}
> {-# LANGUAGE QuasiQuotes #-}

> module Navier.Boundary where

> import Data.Bits
> import Data.Array.Repa

> import Navier.Datadefs
> import Navier.Params

> import Text.Printf
> import Debug.Trace

> apply_boundary_conditions :: Array DIM2 Double -> Array DIM2 Double ->
>                             Array DIM2 Int -> Int -> Int ->
>                             Double -> Double ->
>                             (Array DIM2 Double, Array DIM2 Double)
> apply_boundary_conditions u v flag imax jmax ui vi =
>     let ui' = force $ traverse u id update
>              where update get c@(sh :. j :. i) 
>                     -- Fluid freely flows in from the west
>                     | (i==0) = get (sh :. j :. 1)
>                     -- Fluid freely flows out to the east
>                     | (i==imax) = get (sh :. j :. (imax-1))
>                     | otherwise = get c
>         u' = force $ traverse ui' id update
>              where update get c@(sh :. j :. i) 
>                     -- fluid flows freely in the horizontal direction
>                     | (j==(jmax+1)) = get (sh :. jmax :. i)
>                     | (j==0) = get (sh :. 1 :. i)
>                     | otherwise = get c
>         vi' = force $ traverse v id update
>              where update get c@(sh :. j :. i) 
>                    -- Fluid freely flows in from the west
>                     | (i==0) = get (sh :. j :. 1)
>                    -- Fluid freely flows out to the east
>                     | (i==(imax+1)) = get (sh :. j :. imax)
>                     | otherwise = get c
>         v' = force $ traverse vi' id update
>              where update get c@(sh :. j :. i) 
>                    -- vertical velocity approaches 0 at the north and south
>                     | (j==jmax) = 0.0
>                     | (j==0) = 0.0
>                     | otherwise = get c
>         -- Apply no-slip boundary conditions to cells that are adjacent to
>         -- internal obstacle cells. This forces the u and v velocity to
>         -- tend towards zero in these cells.
>         u'' = force $ traverse2 u' flag (\x y -> x) update
>               where
>                 update getU getFlag c@(sh :. j :. i) =
>                     if (inBounds i j && (getFlag c .&. _bnsew /= 0)) then
>                         case (getFlag c) of
>                           [p|_bn|]  -> -getU (sh :. (j+1) :. i)
>                           [p|_be|]  -> 0.0
>                           [p|_bne|] -> 0.0
>                           [p|_bse|] -> 0.0
>                           [p|_bnw|] -> -getU (sh :. (j+1) :. i)
>                           [p|_bs|]  -> -getU (sh :. (j-1) :. i)
>                           [p|_bsw|] -> -getU (sh :. (j-1) :. i)
>                           _         -> getU c
>                     else
>                         getU c
>         u''' = force $ traverse2 u'' flag (\x y -> x) update
>                where
>                  update getU getFlag c@(sh :. j :. i) =
>                      if (((inBounds i j || i==0) && (i<=(imax-1))) &&
>                          (getFlag (sh :. j :. (i+1)) .&. _bnsew /= 0)) then
>                          case (getFlag (sh :. j :. (i+1))) of
>                            [p|_bn|]  -> -getU (sh :. (j+1) :. i)
>                            [p|_bw|]  -> 0.0
>                            [p|_bne|] -> -getU (sh :. (j+1) :. i)
>                            [p|_bsw|] -> 0.0
>                            [p|_bnw|] -> 0.0
>                            [p|_bs|]  -> -getU (sh :. (j-1) :. i)
>                            [p|_bse|] -> -getU (sh :. (j-1) :. i)
>                            _         -> getU c
>                     else
>                         getU c
>         v'' = force $ traverse2 v' flag (\x y -> x) update
>               where
>                 update getV getFlag c@(sh :. j :. i) =                                                       
>                    if (inBounds i j && (getFlag c .&. _bnsew /= 0)) then
>                      case (getFlag c) of
>                        [p|_bn|] -> 0.0
>                        [p|_be|] -> -getV (sh :. j :. (i+1))
>                        [p|_bne|] -> 0.0
>                        [p|_bse|] -> -getV (sh :. j :. (i+1))
>                        [p|_bnw|] -> 0.0
>                        [p|_bw|] -> -getV (sh :. j :. (i-1))
>                        [p|_bsw|] -> -getV (sh :. j :. (i-1))
>                        _         -> getV c
>                     else
>                         getV c
>         v''' = force $ traverse2 v'' flag (\x y -> x) update
>                where
>                  update getV getFlag c@(sh :. j :. i) =
>                      if (((inBounds i j || j==0) && j<=(jmax-1)) &&
>                          (getFlag (sh :. (j+1) :. i) .&. _bnsew /= 0)) then
>                        case (getFlag (sh :. (j+1) :. i)) of
>                          [p|_be|] -> -getV (sh :. j :. (i+1))
>                          [p|_bs|] -> 0.0
>                          [p|_bne|] -> -getV (sh :. j :. (i+1))
>                          [p|_bse|] -> 0.0
>                          [p|_bsw|] -> 0.0
>                          [p|_bw|] -> -getV (sh :. j :. (i-1))
>                          [p|_bnw|] -> -getV (sh :. j :. (i-1))
>                          _        -> getV c
>                     else
>                         getV c
>         -- Finally, fix the horizontal velocity at the western edge to have
>         -- a continual flow of fluid into the simulation.
>         u4 = force $ traverse u''' id update
>              where update get c@(sh :. j :. i) 
>                     | (i==0 && j>=1 && j<=jmax) = ui
>                     | otherwise = get c
>         v4 = force $ traverse v''' id update
>              where update get c@(sh :. j :. i) 
>                     | (i==0 && j<=jmax) = 2*vi-(get (sh :. j :. 1))
>                     | otherwise = get c
>     in (u4, v4)