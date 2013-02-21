> module Navier.Datadefs where

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