{-# LANGUAGE NoMonomorphismRestriction #-}

module Navier.Init where

import Prelude hiding ( traverse )

import Data.Bits
import Data.Array.Repa

import Navier.Datadefs
import Navier.Params

init_flag :: Int -> Int -> Double -> Double -> (Array U DIM2 Int, Int)
init_flag imax jmax delx dely =
    let
        dims = Z :. (jmax+2) :. (imax+2)
        flag = suspendedComputeP $ fromList dims (Prelude.replicate ((imax+2)*(jmax+2)) (0::Int))
        -- Make a circular obstacle
        mx = 20.0/41.0*(intCast jmax)*dely
        my = mx
        rad1 = 5.0/41.0*(intCast jmax)*dely
        flag' = suspendedComputeP $ traverse flag id update
            where
              update get c@(sh :. j :. i) =
                  if (inBounds i j) then
                            let x = ((intCast i)-0.5)*delx - mx
                                y = ((intCast j)-0.5)*dely - my
                            in
                              if (x*x + y*y <= rad1*rad1) then _cb else _cf
                  else
                      get c
        -- Make the north, sourth, east, and west boundary cells
        flagi' = suspendedComputeP $ traverse flag' id update
            where
              update get c@(sh :. j :. i) =
                  if (j==0) then _cb
                  else if (j==(jmax+1)) then _cb
                       else get c
        flag'' = suspendedComputeP $ traverse flagi' id update
            where
              update get c@(sh :. j :. i) =
                  if (j>=1 && j<=jmax) then
                    if (i==0) then _cb
                    else if (i==(imax+1)) then _cb
                         else get c
                  else
                      get c
        -- Flag boundaries near obstacles
        flag''' = suspendedComputeP $ traverse flag'' id update
             where
               update get c@(sh :. j :. i) =
                   if (inBounds i j && (get c .&. _cf) == 0) then
                       let
                           l = if (get (sh :. j :. (i-1)) .&. _cf /= 0)
                               then _bw else 0
                           r = if (get (sh :. j :. (i+1)) .&. _cf /= 0)
                               then _be else 0
                           b = if (get (sh :. (j-1) :. i) .&. _cf /= 0)
                               then _bs else 0
                           t = if (get (sh :. (j+1) :. i) .&. _cf /= 0)
                               then _bn else 0
                       in
                         (get c) .|. l .|. r .|. b .|. t
                    else
                        get (sh :. j :. i)
        -- Count boundary cells
        ibound = toScalar (foldS total 0 (foldS count 0 (ignoreBoundary flag''')))
                   where count y x = if (x .&. _cf == 0)
                                     then 1+y else y
                         total y x = y + x
   in (flag''', ibound)
