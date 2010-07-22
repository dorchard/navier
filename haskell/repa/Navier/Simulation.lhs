 {-# LANGUAGE QuasiQuotes #-}
 {-# LANGUAGE TypeOperators #-}

> module Navier.Simulation where

> import Navier.Datadefs
> import Navier.Params
> import Data.Array.Repa

> import Debug.Trace
> import Data.Bits

> import Text.Printf

> import Data.Array.Parallel.Base ((:*:)(..))

> set_timestamp_interval del_t (u@Manifest{}) (v@Manifest{}) = 
>     -- del_t satisfying CFL conditions
>     -- else no time stepsize control
>     if (tau >= 1.0e-10) then
>         let
>             umax = foldAll (\x y -> max (abs x) (abs y)) (1.0e-10)
>                 (backpermute (Z :. (jmax-1) :. imax) 
>                              (\(sh :. j :. i) -> (sh :. (j+1) :. i)) u)
>             vmax = foldAll (\x y -> max (abs x) (abs y)) (1.0e-10)
>                 (backpermute (Z :. jmax :. (imax-1))
>                              (\(sh :. j :. i) -> (sh :. j :. (i+1))) v)
>             deltu = delx/umax
>             deltv = dely/vmax
>             deltre = 1.0/(1.0/(delx*delx)+1/(dely*dely))*re/2.0
>         in
>           tau * (if (deltu<deltv) then min deltu deltre
>                                   else min deltv deltre)
>     else
>         del_t

> compute_tentative_velocity  (uA@Manifest{}) (vA@Manifest{})
>                             (fA@Manifest{}) (gA@Manifest{})
>                             (flagA@Manifest{}) del_t =
>    let f'@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) flagA
>             (Data.Array.Repa.zipWith (:*:) uA 
>              (Data.Array.Repa.zipWith (:*:) vA fA))) id updatef
>        flag = fsta
>        u = fsta . snda
>        v = fsta . snda . snda
>        comp = snda . snda. snda
>        updatef get c@(sh :. j :. i) =

>             if ((inBounds i j) && (i<=(imax-1))) then
>              if (((flag $ get c) .&. _cf /= 0)
>                 && ((flag $ get (sh :. j :. (i+1))) .&. _cf /= 0)) then
>                 let

>                     du2dx = ((((u $ get (sh :. j :. i))
>                              +(u $ get (sh :. j :. (i+1))))**2)
>                             +gamma*(abs ((u $ get (sh :. j :. i))
>                                         +(u $ get (sh :. j :. (i+1)))))
>                                   *((u $ get (sh :. j :. i))
>                                    -(u $ get (sh :. j :. (i+1))))
>                             -(((u $ get (sh :. j :. (i-1)))
>                               +(u $ get (sh :. j :. i)))**2)
>                             -gamma*(abs ((u $ get (sh :. j :. (i-1)))
>                                         +(u $ get (sh :. j :. i))))
>                                   *((u $ get (sh :. j :. (i-1)))
>                                    -(u $ get (sh :. j :. i))))/(4.0*delx)
>                     duvdy = ((((v $ get (sh :. j :. i))
>                               +(v $ get (sh :. j :. (i+1))))
>                              *((u $ get (sh :. j :. i))
>                               +(u $ get (sh :. (j+1) :. i))))
>                             +(gamma*(abs ((v $ get (sh :. j :. i))
>                                          +(v $ get (sh :. j :. (i+1)))))
>                                    *((u $ get (sh :. j :. i))
>                                     -(u $ get (sh :. (j+1) :. i))))
>                              -(((v $ get (sh :. (j-1) :. i))
>                                +(v $ get (sh :. (j-1) :. (i+1))))
>                               *((u $ get (sh :. (j-1) :. i))
>                                +(u $ get (sh :. j :. i))))
>                              -(gamma*(abs ((v $ get (sh :. (j-1) :. i))
>                                           +(v $ get (sh :. (j-1) :. (i+1)))))
>                                     *((u $ get (sh :. (j-1) :. i))
>                                      -(u $ get (sh :. j :. i)))))/(4.0*dely)
>                     laplu = ((u $ get (sh :. j :. (i+1)))
>                              -2.0*(u $ get (sh :. j :. i))
>                              +(u $ get (sh :. j :. (i-1))))/delx/delx
>                             +((u $ get (sh :. (j+1) :. i))
>                              -2.0*(u $ get (sh :. j :. i))
>                               + (u $ get (sh :. (j-1) :. i)))/dely/dely
>                 in
>                   (u $ get (sh :. j :. i)) + del_t*(laplu/re-du2dx-duvdy)
>             else
>                 u $ get (sh :. j :. i)
>          else
>              comp $ get (sh :. j :. i)
>        g'@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) flagA
>                       (Data.Array.Repa.zipWith (:*:) uA 
>                        (Data.Array.Repa.zipWith (:*:) vA gA))) id updateg
>        updateg get c@(sh :. j :. i) =
>               if (inBounds i j) && (j<=(jmax-1)) then
>                 if (((flag $ get c) .&. _cf /= 0)
>                   && ((flag $ get (sh :. (j+1) :. i)) .&. _cf /= 0)) then
>                 let
>                     duvdx = ((((u $ get c)+(u $ get (sh :. (j+1) :. i)))
>                              *((v $ get c)+(v $ get (sh :. j :. (i+1)))))
>                             +(gamma*(abs ((u $ get c)+(u $ get (sh :. (j+1) :. i))))
>                                    *((v $ get c)-(v $ get (sh :. j :. (i+1)))))
>                             -(((u $ get (sh :. j :. (i-1)))
>                               +(u $ get (sh :. (j+1) :. (i-1))))
>                               *((v $ get (sh :. j :. (i-1)))+(v $ get c)))
>                             -(gamma*(abs ((u $ get (sh :. j :. (i-1)))
>                                           +(u $ get (sh :. (j+1) :. (i-1)))))
>                                     *((v $ get (sh :. j :. (i-1)))-(v $ get c))))
>                               /(4.0*delx)
>                     dv2dy = ((((v $ get c)+(v $ get (sh :. (j+1) :. i)))**2)
>                             +gamma*(abs ((v $ get c)+(v $ get (sh :. (j+1) :. i))))
>                                   *((v $ get c)-(v $ get (sh :. (j+1) :. i)))
>                             -(((v $ get (sh :. (j-1) :. i))+(v $ get c))**2)
>                             -gamma*(abs ((v $ get (sh :. (j-1) :. i))+(v $ get c)))
>                                   *((v $ get (sh :. (j-1) :. i))-(v $ get c)))/(4.0*dely)
>                     laplv = ((v $ get (sh :. j :. (i+1)))
>                              -2*(v $ get c)
>                              +(v $ get (sh :. j :. (i-1))))/delx/delx +
>                             ((v $ get (sh :. (j+1) :. i))
>                              -2*(v $ get c)
>                              +(v $ get (sh :. (j-1) :. i)))/dely/dely
>                 in
>                   (v $ get c) + del_t*(laplv/re-duvdx-dv2dy)
>                else
>                  v $ get c
>             else
>               comp $ get c

>        f''@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) uA f') id update
>              where
>                update get c@(sh :. j :. i) =
>                    if (j>=1 && j<=jmax) then
>                       if (i==0) then fsta $ get (sh :. j :. 0)
>                       else if (i==imax) then fsta $ get (sh :. j :. imax)
>                       else snda $ get (sh :. j :. i)
>                    else
>                        snda $ get (sh :. j :. i)
>        g''@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) vA g') id update
>              where
>                update get c@(sh :. j :. i) =
>                    if (i>=1 && i<=imax) then
>                       if (j==0) then fsta $ get (sh :. 0 :. i)
>                       else if (j==jmax) then fsta $ get (sh :. jmax :. i)
>                       else snda $ get (sh :. j :. i)
>                    else
>                        snda $ get (sh :. j :. i)
>     in
>       (f'', g'')
>        

> compute_rhs (f@Manifest{}) (g@Manifest{}) (rhs@Manifest{}) 
>              (flag@Manifest{}) del_t = 
>     force $ traverse (Data.Array.Repa.zipWith (:*:) flag
>               (Data.Array.Repa.zipWith (:*:) rhs
>                (Data.Array.Repa.zipWith (:*:) f g))) id update
>     where
>       update get c@(sh :. j :. i) =
>           if (inBounds i j && (((flag $ get c) .&. _cf) /= 0)) then
>              (((f $ get (sh :. j :. i))-(f $ get (sh :. j :. (i-1))))/delx +
>              ((g $ get (sh :. j :. i))-(g $ get (sh :. (j-1) :. i)))/dely)/del_t
>           else
>               rhs $ get (sh :. j :. i)
>             where
>               flag = fsta
>               rhs = fsta . snda
>               f = fsta . snda . snda
>               g = snda . snda . snda


> poisson :: Array DIM2 Double -> Array DIM2 Double -> Array DIM2 Int -> Int ->
>            Double -> (Array DIM2 Double, Double, Int)
> poisson (p@Manifest{}) (rhs@Manifest{}) (flag@Manifest{}) ifluid res = 
>     let
>         rdx2 = 1.0/(delx*delx)
>         rdy2 = 1.0/(dely*dely)
>         beta_2 = -omega/(2.0*(rdx2+rdy2))

>         -- Calculate sum of squares
>         sumSquares = snda $ toScalar (fold total (0 :*: 0.0)
>                                (fold square (0 :*: 0.0) (force $ (ignoreBoundary 
>                        (Data.Array.Repa.zipWith (:*:) flag p)))))
>         total (_ :*: y) (_ :*: x) = (0 :*: (y + x))
>         square (_ :*: y) x  = if ((fsta x) .&. _cf /= 0) then
>                           0 :*: (y+((snda x)**2))
>                       else
>                           0 :*: y

>         p0 = sqrt (sumSquares/((fromIntegral ifluid)::Double))
>         p0' = if (p0 < 0.0001) then 1.0 else p0

>         -- Partial computation of residual
>         computeRes (p@Manifest{}) p0 =
>           let
>             res'@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) flag 
>                                      (Data.Array.Repa.zipWith (:*:) p rhs)) id update
>             res = sumAll $ res'
>             update get c@(sh :. j :. i) =
>                 let
>                    flag = fsta
>                    p = fsta . snda
>                    rhs = snda . snda  
>                 in
>                    if ((inBounds i j) && ((flag $ get c) .&. _cf /= 0)) then
>                      ((((obstacle flag get (sh :. j :. (i+1)))*
>                           ((p $ get (sh :. j :. (i+1)))-(p $ get c)))
>                      -((obstacle flag get (sh :. j :. (i-1)))*
>                           ((p $ get c)-(p $ get (sh :. j :. (i-1))))))*rdx2 
>                     +(((obstacle flag get (sh :. (j+1) :. i))*
>                           ((p $ get (sh :. (j+1) :. i))-(p $ get c)))
>                      -((obstacle flag get (sh :. (j-1) :. i))*
>                           ((p $ get c) - (p $ get (sh :. (j-1) :. i)))))*rdy2
>                     - (rhs $ get c))**2
>                     else
>                         0.0
>           in
>             (sqrt (res/((fromIntegral ifluid)::Double)))/p0

>         -- Red/Black SOR-iteration (compute new p)
>         iterop rb (p@Manifest{}) (rhs@Manifest{}) (flag@Manifest{}) =
>                     force $ traverse (Data.Array.Repa.zipWith (:*:) flag
>                             (Data.Array.Repa.zipWith (:*:) p rhs))
>                             id (update rb)
>         update rb get c@(sh :. j :. i) =
>             let
>               flag = fsta
>               p = fsta . snda
>               rhs = snda . snda
>             in
>               if ((inBounds i j) && (((i+j) `mod` 2) == rb)) then
>                 if ((flag $ get c) == (_cf .|. _bnsew)) then
>                     -- five point star for interior fluid cells
>                    (1.0-omega)*(p $ get c) - 
>                             (beta_2 * 
>                               (((p $ get (sh :. j :. (i+1)))
>                                +(p $ get (sh :. j :. (i-1))))*rdx2
>                              + ((p $ get (sh :. (j+1) :. i))
>                                +(p $ get (sh :. (j-1) :. i)))*rdy2
>                              - (rhs $ get (sh :. j :. i))))
>                 else
>                  -- modified star near boundary
>                  if ((flag $ get c) .&. _cf /= 0) then 
>
>                    let
>                       beta_mod = -omega/((obstacle flag get (sh :. j :. (i+1))+
>                                          (obstacle flag get (sh :. j :. (i-1))))*rdx2
>                                        +((obstacle flag get (sh :. (j+1) :. i))+
>                                          (obstacle flag get (sh :. (j-1) :. i)))*rdy2)
>                    in
>                      (1.0-omega)*(p $ get c) - (beta_mod*(
>                     (((obstacle flag get (sh :. j :. (i+1))*(p $ get (sh :. j :. (i+1)))))
>                     +((obstacle flag get (sh :. j :. (i-1))*(p $ get (sh :. j :. (i-1))))))*rdx2
>                    +(((obstacle flag get (sh :. (j+1) :. i)*(p $ get (sh :. (j+1) :. i))))
>                     +((obstacle flag get (sh :. (j-1) :. i)*(p $ get (sh :. (j-1) :. i)))))*rdy2
>                    - (rhs $ get (sh :. j :. i))))                    
>                  else
>                    p $ get c
>                 else
>                   p $ get c

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

>         -- Iterations
>         sorIterate iter (p@Manifest{}) res = 
>             let
>                 -- red
>                 (p'@Manifest{}) = iterop 0 p rhs flag
>                 -- black
>                 (p''@Manifest{}) = iterop 1 p' rhs flag
>                 -- res
>                 res' = computeRes p'' p0'
>             in
>               if (iter>=itermax) then
>                   (p, res, iter)
>               else 
>                   if (res'<eps) then
>                       (p'', res', iter)
>                   else
>                       sorIterate (iter+1) p'' res'
>                        
>     in
>       sorIterate 0 p res


> update_velocity (u@Manifest{}) (v@Manifest{})
>                 (f@Manifest{}) (g@Manifest{})
>                 (p@Manifest{}) (flag@Manifest{}) del_t =
>     let
>         u'@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) flag
>                               (Data.Array.Repa.zipWith (:*:) p 
>                               (Data.Array.Repa.zipWith (:*:) f u))) id update
>              where
>                update get c@(sh :. j :. i) =
>                    let
>                        flag = fsta
>                        p = fsta . snda
>                        f = fsta . snda . snda
>                        u = snda . snda . snda
>                    in
>                      if (inBounds i j && (i<=(imax-1))
>                                       && ((flag $ get (sh :. j :. i)) .&. _cf /= 0)
>                                       && ((flag $ get (sh :. j :. (i+1))) .&. _cf /= 0)) then
>                          (f $ get c)-((p $ get (sh :. j :. (i+1)))
>                                      -(p $ get (sh :. j :. i)))*del_t/delx
>                      else
>                          u $ get c
>         v'@Manifest{} = force $ traverse (Data.Array.Repa.zipWith (:*:) flag
>                               (Data.Array.Repa.zipWith (:*:) p 
>                               (Data.Array.Repa.zipWith (:*:) g v))) id update
>              where
>                update get c@(sh :. j :. i) =
>                    let
>                        flag = fsta
>                        p = fsta . snda
>                        g = fsta . snda . snda
>                        v = snda . snda . snda
>                    in
>                      if (inBounds i j && (j<=(jmax-1))
>                                       && ((flag $ get (sh :. j :. i)) .&. _cf /= 0)
>                                       && ((flag $ get (sh :. (j+1) :. i)) .&. _cf /= 0)) then
>                          (g $ get c)-((p $ get (sh :. (j+1) :. i))
>                                      -(p $ get (sh :. j :. i)))*del_t/dely
>                      else
>                          v $ get c
>     in
>       (u', v')                 