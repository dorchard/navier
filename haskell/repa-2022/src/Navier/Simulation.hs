{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE TypeOperators #-}

module Navier.Simulation where

import Prelude hiding ( traverse )

import Navier.Datadefs
import Common
import Navier.Params
import qualified Data.Array.Repa as Repa
import Data.Array.Repa
import Data.Array.Repa.Eval ( suspendedComputeP )

import Debug.Trace
import Data.Bits

import Text.Printf

set_timestamp_interval del_t u v =
    -- del_t satisfying CFL conditions
    -- else no time stepsize control
    if (tau >= 1.0e-10) then
        let
            umax = foldAllS (\x y -> max (abs x) (abs y)) (1.0e-10)
                (backpermute (Z :. (jmax-1) :. imax) 
                             (\(sh :. j :. i) -> (sh :. (j+1) :. i)) u)
            vmax = foldAllS (\x y -> max (abs x) (abs y)) (1.0e-10)
                (backpermute (Z :. jmax :. (imax-1))
                             (\(sh :. j :. i) -> (sh :. j :. (i+1))) v)
            deltu = delx/umax
            deltv = dely/vmax
            deltre = 1.0/(1.0/(delx*delx)+1/(dely*dely))*re/2.0
        in
          tau * (if (deltu<deltv) then min deltu deltre
                                  else min deltv deltre)
    else
        del_t

compute_tentative_velocity uA vA fA gA flagA del_t =
   let f' = traverse (Data.Array.Repa.zipWith (,) flagA
            (Data.Array.Repa.zipWith (,) uA
             (Data.Array.Repa.zipWith (,) vA fA))) id updatef
       flag = fst
       u = fst . snd
       v = fst . snd . snd
       comp = snd . snd. snd
       updatef get c@(sh :. j :. i) =

            if ((inBounds i j) && (i<=(imax-1))) then
             if (((flag $ get c) .&. _cf /= 0)
                && ((flag $ get (sh :. j :. (i+1))) .&. _cf /= 0)) then
                let

                    du2dx = ((((u $ get (sh :. j :. i))
                             +(u $ get (sh :. j :. (i+1))))**2)
                            +gamma*(abs ((u $ get (sh :. j :. i))
                                        +(u $ get (sh :. j :. (i+1)))))
                                  *((u $ get (sh :. j :. i))
                                   -(u $ get (sh :. j :. (i+1))))
                            -(((u $ get (sh :. j :. (i-1)))
                              +(u $ get (sh :. j :. i)))**2)
                            -gamma*(abs ((u $ get (sh :. j :. (i-1)))
                                        +(u $ get (sh :. j :. i))))
                                  *((u $ get (sh :. j :. (i-1)))
                                   -(u $ get (sh :. j :. i))))/(4.0*delx)
                    duvdy = ((((v $ get (sh :. j :. i))
                              +(v $ get (sh :. j :. (i+1))))
                             *((u $ get (sh :. j :. i))
                              +(u $ get (sh :. (j+1) :. i))))
                            +(gamma*(abs ((v $ get (sh :. j :. i))
                                         +(v $ get (sh :. j :. (i+1)))))
                                   *((u $ get (sh :. j :. i))
                                    -(u $ get (sh :. (j+1) :. i))))
                             -(((v $ get (sh :. (j-1) :. i))
                               +(v $ get (sh :. (j-1) :. (i+1))))
                              *((u $ get (sh :. (j-1) :. i))
                               +(u $ get (sh :. j :. i))))
                             -(gamma*(abs ((v $ get (sh :. (j-1) :. i))
                                          +(v $ get (sh :. (j-1) :. (i+1)))))
                                    *((u $ get (sh :. (j-1) :. i))
                                     -(u $ get (sh :. j :. i)))))/(4.0*dely)
                    laplu = ((u $ get (sh :. j :. (i+1)))
                             -2.0*(u $ get (sh :. j :. i))
                             +(u $ get (sh :. j :. (i-1))))/delx/delx
                            +((u $ get (sh :. (j+1) :. i))
                             -2.0*(u $ get (sh :. j :. i))
                              + (u $ get (sh :. (j-1) :. i)))/dely/dely
                in
                  (u $ get (sh :. j :. i)) + del_t*(laplu/re-du2dx-duvdy)
            else
                u $ get (sh :. j :. i)
         else
             comp $ get (sh :. j :. i)
       g' = traverse (Data.Array.Repa.zipWith (,) flagA
                      (Data.Array.Repa.zipWith (,) uA
                       (Data.Array.Repa.zipWith (,) vA gA))) id updateg
       updateg get c@(sh :. j :. i) =
              if (inBounds i j) && (j<=(jmax-1)) then
                if (((flag $ get c) .&. _cf /= 0)
                  && ((flag $ get (sh :. (j+1) :. i)) .&. _cf /= 0)) then
                let
                    duvdx = ((((u $ get c)+(u $ get (sh :. (j+1) :. i)))
                             *((v $ get c)+(v $ get (sh :. j :. (i+1)))))
                            +(gamma*(abs ((u $ get c)+(u $ get (sh :. (j+1) :. i))))
                                   *((v $ get c)-(v $ get (sh :. j :. (i+1)))))
                            -(((u $ get (sh :. j :. (i-1)))
                              +(u $ get (sh :. (j+1) :. (i-1))))
                              *((v $ get (sh :. j :. (i-1)))+(v $ get c)))
                            -(gamma*(abs ((u $ get (sh :. j :. (i-1)))
                                          +(u $ get (sh :. (j+1) :. (i-1)))))
                                    *((v $ get (sh :. j :. (i-1)))-(v $ get c))))
                              /(4.0*delx)
                    dv2dy = ((((v $ get c)+(v $ get (sh :. (j+1) :. i)))**2)
                            +gamma*(abs ((v $ get c)+(v $ get (sh :. (j+1) :. i))))
                                  *((v $ get c)-(v $ get (sh :. (j+1) :. i)))
                            -(((v $ get (sh :. (j-1) :. i))+(v $ get c))**2)
                            -gamma*(abs ((v $ get (sh :. (j-1) :. i))+(v $ get c)))
                                  *((v $ get (sh :. (j-1) :. i))-(v $ get c)))/(4.0*dely)
                    laplv = ((v $ get (sh :. j :. (i+1)))
                             -2*(v $ get c)
                             +(v $ get (sh :. j :. (i-1))))/delx/delx +
                            ((v $ get (sh :. (j+1) :. i))
                             -2*(v $ get c)
                             +(v $ get (sh :. (j-1) :. i)))/dely/dely
                in
                  (v $ get c) + del_t*(laplv/re-duvdx-dv2dy)
               else
                 v $ get c
            else
              comp $ get c

       f'' = suspendedComputeP $ traverse (Data.Array.Repa.zipWith (,) uA f') id update
             where
               update get c@(sh :. j :. i) =
                   if (j>=1 && j<=jmax) then
                      if (i==0) then fst $ get (sh :. j :. 0)
                      else if (i==imax) then fst $ get (sh :. j :. imax)
                      else snd $ get (sh :. j :. i)
                   else
                       snd $ get (sh :. j :. i)
       g'' = suspendedComputeP $ traverse (Data.Array.Repa.zipWith (,) vA g') id update
             where
               update get c@(sh :. j :. i) =
                   if (i>=1 && i<=imax) then
                      if (j==0) then fst $ get (sh :. 0 :. i)
                      else if (j==jmax) then fst $ get (sh :. jmax :. i)
                      else snd $ get (sh :. j :. i)
                   else
                       snd $ get (sh :. j :. i)
    in
      (f'', g'')

compute_rhs f g rhs flag del_t =
    suspendedComputeP $ traverse (Data.Array.Repa.zipWith (,) flag
              (Data.Array.Repa.zipWith (,) rhs
               (Data.Array.Repa.zipWith (,) f g))) id update
    where
      update get c@(sh :. j :. i) =
          if (inBounds i j && (((flag $ get c) .&. _cf) /= 0)) then
             (((f $ get (sh :. j :. i))-(f $ get (sh :. j :. (i-1))))/delx +
             ((g $ get (sh :. j :. i))-(g $ get (sh :. (j-1) :. i)))/dely)/del_t
          else
              rhs $ get (sh :. j :. i)
            where
              flag = fst
              rhs = fst . snd
              f = fst . snd . snd
              g = snd . snd . snd


poisson
    :: Array U DIM2 Double -> Array U DIM2 Double -> Array U DIM2 Int
    -> Int -> Double
    -> (Array U DIM2 Double, Double, Int)
poisson p rhs flag ifluid res =
    let
        rdx2 = 1.0/(delx*delx)
        rdy2 = 1.0/(dely*dely)
        beta_2 = -omega/(2.0*(rdx2+rdy2))

        -- Calculate sum of squares
        sumSquares :: Double
        sumSquares =
            snd $ toScalar $ sumSquaresCompute sumSquaresInner
        sumSquaresInner :: Array D DIM2 (Int, Double)
        sumSquaresInner = ignoreBoundary $ Data.Array.Repa.zipWith (,) flag p
        sumSquaresCompute :: Array D DIM2 (Int, Double) -> Array U DIM0 (Int, Double)
        sumSquaresCompute = foldS total (0, 0.0) . foldS square (0, 0.0)

        total (_, y) (_, x) = (0, (y + x))
        square (_, y) x  =
            if   ((fst x) .&. _cf /= 0)
            then (0, (y+((snd x)**2)))
            else (0, y)

        p0 = sqrt (sumSquares/((fromIntegral ifluid)::Double))
        p0' = if (p0 < 0.0001) then 1.0 else p0

        -- Partial computation of residual
        computeRes p p0 =
          let
            res' = traverse (Data.Array.Repa.zipWith (,) flag
                                     (Data.Array.Repa.zipWith (,) p rhs)) id update
            res = sumAllS $ res'
            update get c@(sh :. j :. i) =
                let
                   flag = fst
                   p = fst . snd
                   rhs = snd . snd
                in
                   if ((inBounds i j) && ((flag $ get c) .&. _cf /= 0)) then
                     ((((obstacle flag get (sh :. j :. (i+1)))*
                          ((p $ get (sh :. j :. (i+1)))-(p $ get c)))
                     -((obstacle flag get (sh :. j :. (i-1)))*
                          ((p $ get c)-(p $ get (sh :. j :. (i-1))))))*rdx2 
                    +(((obstacle flag get (sh :. (j+1) :. i))*
                          ((p $ get (sh :. (j+1) :. i))-(p $ get c)))
                     -((obstacle flag get (sh :. (j-1) :. i))*
                          ((p $ get c) - (p $ get (sh :. (j-1) :. i)))))*rdy2
                    - (rhs $ get c))**2
                    else
                        0.0
          in
            (sqrt (res/((fromIntegral ifluid)::Double)))/p0

        -- Red/Black SOR-iteration (compute new p)
        iterop rb p rhs flag =
                    traverse (Data.Array.Repa.zipWith (,) flag
                            (Data.Array.Repa.zipWith (,) p rhs))
                            id (update rb)
        update rb get c@(sh :. j :. i) =
            let
              flag = fst
              p = fst . snd
              rhs = snd . snd
              -- CSE the central element access
              !_c = get $ c
            in
              if ((inBounds i j) && (((i+j) `mod` 2) == rb)) then
                if ((flag _c) == (_cf .|. _bnsew)) then
                    -- five point star for interior fluid cells
                   (1.0-omega)*(p _c) - 
                            (beta_2 * 
                              (((p $ get (sh :. j :. (i+1)))
                               +(p $ get (sh :. j :. (i-1))))*rdx2
                             + ((p $ get (sh :. (j+1) :. i))
                               +(p $ get (sh :. (j-1) :. i)))*rdy2
                             - (rhs _c)))
                else
                 -- modified star near boundary
                 if ((flag _c) .&. _cf /= 0) then 
                   let
                      beta_mod = -omega/((obstacle flag get (sh :. j :. (i+1))+
                                         (obstacle flag get (sh :. j :. (i-1))))*rdx2
                                       +((obstacle flag get (sh :. (j+1) :. i))+
                                         (obstacle flag get (sh :. (j-1) :. i)))*rdy2)
                   in
                     (1.0-omega)*(p _c) - (beta_mod*(
                    (((obstacle flag get (sh :. j :. (i+1))*(p $ get (sh :. j :. (i+1)))))
                    +((obstacle flag get (sh :. j :. (i-1))*(p $ get (sh :. j :. (i-1))))))*rdx2
                   +(((obstacle flag get (sh :. (j+1) :. i)*(p $ get (sh :. (j+1) :. i))))
                    +((obstacle flag get (sh :. (j-1) :. i)*(p $ get (sh :. (j-1) :. i)))))*rdy2
                   - (rhs _c)))                    
                 else
                   p _c
                else
                  p _c

{-

                ----- CSE version (actually slower) 

                 let
                   _l = get (sh :. j :. (i-1)) 
                   _r = get (sh :. j :. (i+1))
                   _t = get (sh :. (j-1) :. i)
                   _b = get (sh :. (j+1) :. i)
                   beta_mod = -omega/(((obs $ flag _r)+
                                         (obs $ flag _l))*rdx2
                                       +((obs $ flag _b)+
                                         (obs $ flag _t))*rdy2)
                  in
                    (1.0-omega)*(p $ get c) - (beta_mod*(
                   (
                    (((obs $ flag _r)*p _r)
                    +((obs $ flag _l)*p _l))*rdx2
                   +(((obs $ flag _b)*p _b)
                    +((obs $ flag _t)*p _t))*rdy2
                   -(rhs $ get c))))
                else
                  p $ get c
               else
                 p $ get c

-}

        -- Iterations
        sorIterate :: Int -> Array U DIM2 Double -> Double -> (Array U DIM2 Double, Double, Int)
        sorIterate iter p res =
            let
                -- red
                p' = iterop 0 p rhs flag
                -- black
                p'' = iterop 1 p' rhs flag
                -- res
                res' = computeRes p'' p0'
                p''' = suspendedComputeP p''
            in
              if (iter>=itermax) then
                  (p, res, iter)
              else 
                  if (res'<eps) then
                      (p''', res', iter)
                  else
                      sorIterate (iter+1) p''' res'

    in
      sorIterate 0 p res


update_velocity u v f g p flag del_t =
    let
        u' = suspendedComputeP $ traverse (Data.Array.Repa.zipWith (,) flag
                              (Data.Array.Repa.zipWith (,) p
                              (Data.Array.Repa.zipWith (,) f u))) id update
             where
               update get c@(sh :. j :. i) =
                   let
                       flag = fst
                       p = fst . snd
                       f = fst . snd . snd
                       u = snd . snd . snd
                   in
                     if (inBounds i j && (i<=(imax-1))
                                      && ((flag $ get (sh :. j :. i)) .&. _cf /= 0)
                                      && ((flag $ get (sh :. j :. (i+1))) .&. _cf /= 0)) then
                         (f $ get c)-((p $ get (sh :. j :. (i+1)))
                                     -(p $ get (sh :. j :. i)))*del_t/delx
                     else
                         u $ get c
        v' = suspendedComputeP $ traverse (Data.Array.Repa.zipWith (,) flag
                              (Data.Array.Repa.zipWith (,) p
                              (Data.Array.Repa.zipWith (,) g v))) id update
             where
               update get c@(sh :. j :. i) =
                   let
                       flag = fst
                       p = fst . snd
                       g = fst . snd . snd
                       v = snd . snd . snd
                   in
                     if (inBounds i j && (j<=(jmax-1))
                                      && ((flag $ get (sh :. j :. i)) .&. _cf /= 0)
                                      && ((flag $ get (sh :. (j+1) :. i)) .&. _cf /= 0)) then
                         (g $ get c)-((p $ get (sh :. (j+1) :. i))
                                     -(p $ get (sh :. j :. i)))*del_t/dely
                     else
                         v $ get c
    in
      (u', v')