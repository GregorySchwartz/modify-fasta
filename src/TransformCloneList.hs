-- TransformCloneList module.
-- By Gregory W. Schwartz
--
-- Collection of functions that transform a clone entry in some way

{-# LANGUAGE BangPatterns, OverloadedStrings #-}

module TransformCloneList ( onlyMutations
                          , frequentMutations
                          ) where

-- Built-in
import Data.Maybe
import Data.List
import qualified Data.Text as T
import qualified Data.Map.Strict as Map
import Data.Tuple
import Data.Monoid
import Control.Arrow (second)

-- Cabal
import Data.Fasta.Text
import qualified Data.List.Split as Split

-- Local
import Types
import Utility
import Diversity

-- | Return True if there are no gap characters in the text
noGaps :: T.Text -> Bool
noGaps = not . any (\x -> x == '.' || x == '-') . T.unpack

-- Replace codons that have more than CodonMut mutations (make them "---"
-- codons) and don't have gaps in them.
onlyMutations :: CodonMut -> T.Text -> T.Text -> CloneEntry -> CloneEntry
onlyMutations codonMut codonMutType mutType = newEntry
  where
    newEntry (!germ, !fseqs)  = (germ, map (removeCodon germ) fseqs)
    removeCodon germ fseq     = fseq { fastaSeq = remove (fastaSeq germ)
                                                . fastaSeq
                                                $ fseq
                                     }
    remove germSeq            = mconcat
                              . map (snd . replaceCodon)
                              . zip (codonSplit germSeq)
                              . codonSplit
    replaceCodon (x, y)
        | (codonMutOp codonMutType) (hamming x y) codonMut
       && isMutType (T.toUpper mutType) x y
       && noGaps x
       && noGaps y  = (x, y)
        | otherwise = ("---", "---")
    codonSplit                   = fullCodon . T.chunksOf 3
    fullCodon                    = filter ((== 3) . T.length)
    codonMutOp ">" = (>)
    codonMutOp "<" = (<)
    codonMutOp "=" = (==)
    isMutType "REPLACEMENT" x y = codon2aa x /= codon2aa y
    isMutType "SILENT" x y      = codon2aa x == codon2aa y
    isMutType "ALL" _ _         = True
    isMutType _ _ _             = error "Unknown mutation type"

-- Only include codons containing mutations found in a certain number of
-- mutants
frequentMutations :: Maybe Int
                  -> Maybe Int
                  -> Maybe Double
                  -> CloneEntry
                  -> CloneEntry
frequentMutations minSeqs mutCount mutPercent entry@(!germline, !fseqs) =
    newEntry
  where
    newEntry         = (germline, map removeCodon fseqs)
    removeCodon fseq = fseq { fastaSeq = replace fseq }
    replace          = rejoinCodons
                     . replaceCodon minSeqs mutCount mutPercent numSeqs countMap
                     . positionalCodons germline
    rejoinCodons     = T.pack . concatMap (map (snd . snd))
    countMap         = getCountMap germline fseqs
    numSeqs          = length fseqs

-- | Replace codons that are not valid with a gap "---". Important to note
-- that the predicates return False if we have what we want because the
-- monoid instance of Any has True trumping False, but we want False to
-- trump True, so invert it all and invert it back at the end.
replaceCodon :: Maybe Int
             -> Maybe Int
             -> Maybe Double
             -> Int
             -> CountMap
             -> CodonMutations
             -> CodonMutations
replaceCodon minSeqs mutCount mutPercent numSeqs countMap = map replace
  where
    replace x = if (any isValid . filter isMutation $ x)
                && (not . any (isGapMut . snd) $ x)
                    then x
                    else map ((second . second) (const '-')) x
    isGapMut (x, y) = x == '.' || x == '-' || y == '.' || y == '-'
    isValid x = case Map.lookup x countMap of
            (Just count) -> not
                          . fromMaybe False
                          . fmap getAny
                          $ minSeqsValid
                         <> mutCountValid count
                         <> mutPercentValid count
            Nothing -> error ("Mutation not found: " ++ show x)
    minSeqsValid           = fmap (Any . not . (>=) numSeqs) minSeqs
    mutCountValid count    = case mutCount of
                                 (Just 0) -> fmap (Any . not . (==) count)
                                           $ Just numSeqs
                                 (Just x) -> Just . Any . not . (>=) count $ x
                                 Nothing  -> Nothing
    mutPercentValid count  =
        fmap
        (Any . not . (>=) ((fromIntegral count / fromIntegral numSeqs) * 100))
        mutPercent
    isMutation (_, (x, y)) = x /= y

-- | Get the complete mutation codons with a position provided for each
-- nucleotide
positionalCodons :: Germline -> FastaSequence -> CodonMutations
positionalCodons germline = fullCodon . codonSplit . T.unpack . fastaSeq
  where
    codonSplit = fullCodon
               . Split.chunksOf 3
               . zip [1..]
               . zip (T.unpack . fastaSeq $ germline)
    fullCodon  = filter ((== 3) . length)

-- | Get the number of times each mutation appears
getCountMap :: Germline
            -> [FastaSequence]
            -> CountMap
getCountMap (FastaSequence { fastaSeq = germlineSeq }) =
    Map.unionsWith (+) . map initializeMaps
  where
    initializeMaps = Map.fromList
                   . map swap
                   . zip [1,1..]
                   . zip [1..]
                   . T.zip germlineSeq
                   . fastaSeq

