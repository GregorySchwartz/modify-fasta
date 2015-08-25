-- Utility module
-- By Gregory W. Schwartz
--
-- Collects utility functions for the main files

module Utility ( addLengthHeader
               , addFillerGermlines
               , replaceChars
               ) where

-- Built-in
import qualified Data.Map as M

-- Cabal
import Data.Fasta.String

-- | Adds the length of a sequence to the header of that sequence
addLengthHeader :: FastaSequence -> FastaSequence
addLengthHeader fSeq = fSeq { fastaHeader = fastaHeader fSeq
                                         ++ "|"
                                         ++ (show . length . fastaSeq $ fSeq)
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

-- Replace characters in the first string with another in the second string
-- if they are equal to a certain character and they aren't replaced with
-- a gap.
replaceChars :: Char -> String -> String -> String
replaceChars c = zipWithRetain changeChar
  where
    changeChar a b = if a == c && b `notElem` ".-" then b else a
