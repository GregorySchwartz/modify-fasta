-- FilterFastaList module.
-- By Gregory W. Schwartz
--
-- Collection of functions for the filtering of a pipesFasta

{-# LANGUAGE BangPatterns, OverloadedStrings, FlexibleContexts #-}

module FilterCloneList ( filterHighlyMutatedEntry
                       ) where

-- Built in
import Data.List
import Data.Maybe
import Data.Either
import Text.Regex.TDFA
import Text.Regex.TDFA.Text
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text

-- Local
import Types

-- | Remove highly mutated sequences (sequences with more than a third of
-- their sequence being mutated).
filterHighlyMutatedEntry :: GeneticUnit
                         -> CodonTable
                         -> CloneEntry
                         -> CloneEntry
filterHighlyMutatedEntry !genUnit !table = newEntry
  where
    newEntry (!germline, !fseqs) = ( germline
                                   , map snd
                                   . filter (not . fst)
                                   . rights
                                   . assignMutated germline
                                   $ fseqs
                                   )
    assignMutated k              = map (isHighlyMutated k)
    isHighlyMutated !k !x        =
        case (readSeq genUnit k, readSeq genUnit x) of
            ((Right a), (Right b)) -> (\n -> Right (n, b))
                                    $ ( (fromIntegral (T.length (fastaSeq a)) :: Double)
                                      / 3 )
                                   <= ( ( genericLength
                                        . realMutations (fastaSeq a)
                                        $ fastaSeq b ) )
            ((Left a), (Right _)) -> error ("Error in germline: " ++ T.unpack a)
            ((Right _), (Left b)) -> error ("Error in sequence: " ++ T.unpack b)
            ((Left a), (Left b))  -> error (unwords [ "Error in sequence:"
                                                    , T.unpack b
                                                    , "with germline:"
                                                    , T.unpack a ] )
    realMutations k x   = filterCodonMutStab (\(!y, !z) -> y /= z)
                        . map snd
                        . mutation k
                        $ x
    filterCodonMutStab isWhat = filter (filterRules genUnit isWhat)
    filterRules AminoAcid isWhat x = isWhat x
                                  && not (inTuple '-' x)
                                  && not (inTuple '.' x)
                                  && not (inTuple '~' x)
    filterRules Nucleotide isWhat x = isWhat x
                                   && not (inTuple '-' x)
                                   && not (inTuple '.' x)
                                   && not (inTuple '~' x)
                                   && not (inTuple 'N' x)
    inTuple c (x, y)
        | c == x || c == y = True
        | otherwise        = False
    mutation x y        = zip [1..] . T.zip x $ y
    readSeq Nucleotide x = Right x
    readSeq AminoAcid x  = customTranslate table 1 x
