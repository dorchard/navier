module Common where

import Data.Array.Repa

-- TODO likely a strict pair so.
data a :*: b = !a :*: !b
infixr 9 :*:

fsta :: a :*: b -> a
fsta (a :*: _) = a

snda :: a :*: b -> b
snda (_ :*: b) = b

toScalar :: Array U DIM0 Int -> Int
toScalar arr = index arr Z
