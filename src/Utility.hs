-- Utility module
-- By Gregory W. Schwartz
--
-- Collects utility functions for the main files

{-# LANGUAGE OverloadedStrings, ViewPatterns #-}

module Utility ( addLengthHeader
               , addFillerGermlines
               , replaceChars
               ) where

-- Built-in
import qualified Data.Map as M
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text
import TextShow

-- | Adds the length of a sequence to the header of that sequence
addLengthHeader :: FastaSequence -> FastaSequence
addLengthHeader fSeq = fSeq { fastaHeader = fastaHeader fSeq
                                  `mappend` "|"
                                  `mappend` (showt . T.length . fastaSeq $ fSeq)
                            }

-- | Adds filler germlines to normal fasta files
addFillerGermlines :: [FastaSequence] -> CloneMap
addFillerGermlines = M.fromList . labelGermlines . map insertDummy
  where
    labelGermlines  = map (\(x, (y, z)) -> ((x, y), z)) . zip [0..]
    insertDummy x   = (dummy, [x])
    dummy = FastaSequence {fastaHeader = "filler", fastaSeq = "---"}

-- Like zipWith, but if one if one list is longer than the other than use
-- the remaining, needs to be the same type
zipWithRetain :: (a -> a -> a) -> [a] -> [a] -> [a]
zipWithRetain _ [] [] = []
zipWithRetain _ xs [] = xs
zipWithRetain _ [] ys = ys
zipWithRetain f (x:xs) (y:ys) = f x y : zipWithRetain f xs ys

-- Like zipWithRetain, but for text
zipWithRetainText :: (Char -> Char -> Char) -> T.Text -> T.Text -> T.Text
zipWithRetainText _ (T.uncons -> Nothing) (T.uncons -> Nothing) = T.empty
zipWithRetainText _ xs (T.uncons -> Nothing) = xs
zipWithRetainText _ (T.uncons -> Nothing) ys = ys
zipWithRetainText f (T.uncons -> Just (x, xs)) (T.uncons -> Just (y, ys))
    = f x y `T.cons` zipWithRetainText f xs ys

-- Replace characters in the first string with another in the second string
-- if they are equal to a certain character and they aren't replaced with
-- a gap.
replaceChars :: Char -> T.Text -> T.Text -> T.Text
replaceChars c = zipWithRetainText changeChar
  where
    changeChar a b = if a == c && (not . T.isInfixOf (T.singleton b)) ".-"
                        then b
                        else a
