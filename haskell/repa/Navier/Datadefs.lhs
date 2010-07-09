> {-# LANGUAGE NoMonomorphismRestriction #-}
> {-# LANGUAGE TemplateHaskell #-}

> module Navier.Datadefs where

> import Data.Bits
> import Language.Haskell.TH
> import Language.Haskell.TH.Quote

> import Data.Array.Parallel.Base ((:*:)(..))

> import Debug.Trace

> fsta (a :*: b) = a
> snda (a :*: b) = b

> _cb = 0x0000
> _bn = 0x0001
> _bs = 0x0002
> _bw = 0x0004
> _be = 0x0008
> _bnw = _bn .|. _bw
> _bsw = _bs .|. _bw
> _bne = _bn .|. _be
> _bse = _bs .|. _be
> _bnsew = _bn .|. _bs .|. _be .|. _bw
> _cf = 0x0010

> p :: QuasiQuoter
> p = QuasiQuoter undefined constPat undefined undefined

> constPat :: String -> PatQ
> constPat varName = litP $ integerL $ case varName of
>                      "_be" -> _be
>                      "_bn" -> _bn
>                      "_bs" -> _bs
>                      "_bw" -> _bw
>                      "_bne" -> _bne
>                      "_bse" -> _bse
>                      "_bsw" -> _bsw
>                      "_bnw" -> _bnw
>                      _ -> trace "FAIL!" 0
>                               
