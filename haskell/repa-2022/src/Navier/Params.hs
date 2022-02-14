module Navier.Params where

import Prelude hiding ( traverse )

import Data.Bits
import Data.Array.Repa
import Navier.Datadefs
import Data.Vector.Unboxed.Base ( Unbox )

-- boundary flags
xlength = 22.0       -- Width of simulated domain
ylength = 4.1        -- Height of simulated domain
imax = 660           -- Number of cells horizontally
jmax = 120           -- Number of cells vertically

delx = xlength/660.0
dely = ylength/120.0

t_end = 40.0         -- Simultion runtime
del_t = 0.003        -- (default) duration of each timestep
tau = 0.5            -- safety factor for timestep control

itermax = 10::Int   -- max number of SOR iterations
eps = 0.001          -- stopping error threshold for SOR
omega = 1.7          -- relaxation parameter for SOR
gamma = 0.9          -- upwind differencing factor in PDE discretisation

re = 150.0           -- Reynolds number
ui = 1.0             -- Initial X velocity (West-east flow)
vi = 0.0             -- Initial Y velocity

{-# INLINE inBounds #-}
inBounds :: Int -> Int -> Bool
inBounds i j = (i>0) && (i<(imax+1)) && (j>0) && (j<(jmax+1))

intCast :: Int -> Double
intCast = fromInteger . toInteger

ignoreBoundary :: (Unbox a, Source r a) => Array r DIM2 a -> Array D DIM2 a
ignoreBoundary arr = let Z :. height :. width = extent arr
                     in backpermute (Z :. (height-2) :. (width-2))
                           (\(sh :. y :. x) -> (sh :. (y+1) :. (x+1))) arr

{-# INLINE obs #-}
obs x = if (x .&. _cf /= 0) then 1 else 0

{-# INLINE obstacle #-}
obstacle acc get i = if ((acc $ get i) .&. _cf /= 0) then 1 else 0

checksum_char :: (Source r Int, Monad m) => Array r DIM2 Int -> m Int
checksum_char arr = sumAllP $ traverse arr id update
                    where update get c@(sh :. j :. i) =
                              (get c)*(i)

checksum :: (Source r Double, Monad m) => Array r DIM2 Double -> m Double
checksum arr = sumAllP $ traverse arr id update
                    where update get c@(sh :. j :. i) =
                              (get c)*(intCast $ (i*jmax)+j)
