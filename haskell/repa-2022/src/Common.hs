module Common where

import Data.Array.Repa
import Data.Vector.Unboxed

toScalar :: Source r a => Array r DIM0 a -> a
toScalar arr = index arr Z
