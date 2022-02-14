{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE TemplateHaskell #-}

module Navier.Datadefs where

import Prelude hiding ( (++) )

import Data.Bits
import Language.Haskell.TH
import Language.Haskell.TH.Quote
import Data.Array.Repa

import Debug.Trace

_cb = 0x0000
_bn = 0x0001
_bs = 0x0002
_bw = 0x0004
_be = 0x0008
_bnw = _bn .|. _bw
_bsw = _bs .|. _bw
_bne = _bn .|. _be
_bse = _bs .|. _be
_bnsew = _bn .|. _bs .|. _be .|. _bw
_cf = 0x0010

p :: QuasiQuoter
p = QuasiQuoter undefined constPat undefined undefined

constPat :: String -> PatQ
constPat varName = litP $ integerL $ case varName of
                     "_be" -> _be
                     "_bn" -> _bn
                     "_bs" -> _bs
                     "_bw" -> _bw
                     "_bne" -> _bne
                     "_bse" -> _bse
                     "_bsw" -> _bsw
                     "_bnw" -> _bnw
                     _ -> trace "FAIL!" 0

{-

stencil :: QuasiQuoter
stencil = QuasiQuoter stencilExp stencilPat undefined stencilDecl

tie _ [] = []
tie [] _ = []
tie (x:xs) (y:ys) = (x++y):(tie xs ys)

pos :: [(String, Int, Int)]
pos = [("_tl", -1, -1), ("_t", 0, -1), ("_tr", 1, -1),
       ("_l", -1, 0), ("_c", 0, 0), ("_r", 1, 0),
       ("_bl", -1, 1), ("_b", 0, 1), ("_br", 1, 1)]

stencilDecl :: String -> Q [Dec]
stencilDecl varsString =
    do vars <- return $ words varsString
       mapM (binder vars) pos

binder vars (l, x, y) =
       let pat = Prelude.map (varP . mkName)  (tie vars (repeat l))
           exp = appE (varE $ mkName "get")
                  [| ($(varE $ mkName "sh")) :.
                     ($(varE $ mkName "j") + y) :.
                     ($(varE $ mkName "i") + x) |]
       in valD (tupP pat) (normalB exp) []

stencilPat :: String -> Q Pat
stencilPat varsString
             = do vars <- return $ words varsString
                  tupP [tupP (Prelude.map (mkPat) pos),
                        tupP (Prelude.map (varP . mkName) vars),
                        tupP (Prelude.concatMap (mkPats vars) pos)]

mkPat (l, x, y) = varP $ mkName l

mkPats vars (l, _, _) = Prelude.map (varP . mkName) (tie vars (repeat l))

inter [v] l = varP $ mkName (v++l)
inter (v:vs) l = infixP (varP $ mkName (v++l)) (mkName ":*:") (inter vs l)

--Prelude.map (varP . mkName) (tie vars (repeat l))

stencilExp :: String -> Q Exp
stencilExp varsString = do vars <- return $ words varsString
                           tupE [tupE (Prelude.map mkExps pos),
                                 tupE (accessors vars),
                                 tupE (Prelude.concatMap (mappings vars) pos)]

mkExps (_, x, y) =
             appE (varE $ mkName "get")
                  [| ($(varE $ mkName "sh")) :. 
                     ($(varE $ mkName "j") + y) :.
                     ($(varE $ mkName "i") + x) |]

accessors vars = Prelude.map conv (accs vars)

--(iterate add1 [False]))

mappings vars (l, x, y) = Prelude.map (\v -> appE (varE $ mkName v)
                                             (varE $ mkName l)) vars

conv xs = Prelude.foldr1 (\x y -> infixE (Just x) (varE $ mkName ".") (Just y))
             (Prelude.map (varE . mkName) xs)

add1 :: [Bool] -> [Bool]
add1 [] = [True]
add1 (x:xs) = if x then (not x):(add1 xs) else True:xs

accs vars = accs' (tail vars) ["fsta"]
accs' [x] acc = [acc, "snda":(tail acc)]
accs' (x:xs) acc = acc:(accs' xs (acc++["snda"]))

-}
