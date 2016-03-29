-- Utility module
-- By Gregory W. Schwartz
--
-- Collects utility functions for the main files

{-# LANGUAGE BangPatterns, OverloadedStrings, ViewPatterns #-}

module Utility ( addLengthHeader
               , addMutationsHeader
               , addFillerGermlines
               , replaceChars
               , getField
               , fromEither
               ) where

-- Built-in
import qualified Data.Map as M
import Data.Monoid

-- Cabal
import qualified Data.Text as T
import Data.Fasta.Text
import TextShow

-- Local
import Types

-- | Adds the length of a sequence to the header of that sequence
addLengthHeader :: FastaSequence -> FastaSequence
addLengthHeader fSeq = fSeq { fastaHeader = fastaHeader fSeq
                                         <> "|"
                                         <> (showt . T.length . fastaSeq $ fSeq)
                            }

-- | Adds the mutations of a sequence to the header of that sequence
addMutationsHeader :: Bool -> Field -> FastaSequence -> FastaSequence
addMutationsHeader aaFlag field fSeq =
    fSeq { fastaHeader = fastaHeader fSeq
                      <> "|"
                      <> ( printMutations
                         . getMutations (fastaSeq germline)
                         . fastaSeq
                         $ fSeq
                         )
         }
  where
    germline = if aaFlag then fromEither (translate 1 otherSeq) else otherSeq
    otherSeq = FastaSequence {fastaHeader = "", fastaSeq = getField field fSeq}

-- | Print the mutations
printMutations :: [(Position, (Char, Char))] -> T.Text
printMutations = T.intercalate "/"
               . map (\(!p, (!x, !y)) -> showt p <> T.pack ['-', x, '-', y])

-- | Filter for the true mutations
getMutations :: T.Text -> T.Text -> [(Position, (Char, Char))]
getMutations xs = filter (\x -> isDiff x && noGaps x) . getDiff xs
  where
    isDiff (_, (!x, !y)) = x /= y
    noGaps (_, !x)       = not . any (flip inTuple x) $ ("-.~" :: String)

-- | Returns the difference between two texts
getDiff :: T.Text -> T.Text -> [(Position, (Char, Char))]
getDiff xs = zip [1..] . T.zip xs

-- | Sees if an element is in a tuple
inTuple :: (Eq a) => a -> (a, a) -> Bool
inTuple c (!x, !y) = c == x || c == y

-- | Adds filler germlines to normal fasta files
addFillerGermlines :: [FastaSequence] -> CloneMap
addFillerGermlines = M.fromList . labelGermlines . map insertDummy
  where
    labelGermlines  = map (\(x, (y, z)) -> ((x, y), z)) . zip [0..]
    insertDummy x   = (dummy, [x])
    dummy = FastaSequence {fastaHeader = "filler", fastaSeq = "---"}

-- | Like zipWith, but if one if one list is longer than the other than use
-- the remaining, needs to be the same type
zipWithRetain :: (a -> a -> a) -> [a] -> [a] -> [a]
zipWithRetain _ [] [] = []
zipWithRetain _ xs [] = xs
zipWithRetain _ [] ys = ys
zipWithRetain f (x:xs) (y:ys) = f x y : zipWithRetain f xs ys

-- | Like zipWithRetain, but for text
zipWithRetainText :: (Char -> Char -> Char) -> T.Text -> T.Text -> T.Text
zipWithRetainText _ (T.uncons -> Nothing) (T.uncons -> Nothing) = T.empty
zipWithRetainText _ xs (T.uncons -> Nothing) = xs
zipWithRetainText _ (T.uncons -> Nothing) ys = ys
zipWithRetainText f (T.uncons -> Just (x, xs)) (T.uncons -> Just (y, ys))
    = f x y `T.cons` zipWithRetainText f xs ys

-- | Replace characters in the first string with another in the second string
-- if they are equal to a certain character and they aren't replaced with
-- a gap.
replaceChars :: Char -> T.Text -> T.Text -> T.Text
replaceChars c = zipWithRetainText changeChar
  where
    changeChar a b = if a == c && (not . T.isInfixOf (T.singleton b)) ".-"
                        then b
                        else a

-- | Get the field of a fasta sequence, 1 indexed split by "|"
getField :: Int -> FastaSequence -> T.Text
getField f fs = (T.splitOn "|" . fastaHeader $ fs) !! (f - 1)

-- | Error for left
fromEither :: Either T.Text b -> b
fromEither (Right x)     = x
fromEither (Left x)      = error . T.unpack $ x
