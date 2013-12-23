-- FilterCloneMap module.
-- By G.W. Schwartz
--
-- Collection of functions for the filtering of a CloneMap

module FilterCloneMap where

-- Built in
import Data.List
import qualified Data.Map as M

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types
import Diversity
import Translation

-- Replace codons that have more than CodonMut mutations (make them "---"
-- codons).
removeCodonMutCount :: CodonMut -> CloneMap -> CloneMap
removeCodonMutCount codonMut = M.mapWithKey mapRemove
  where
    mapRemove (pos, germ) clones = map (removeCodon germ) clones
    removeCodon germ clone       = clone { fastaSeq
                                         = remove (fastaSeq germ)
                                         . fastaSeq $ clone }
    remove germSeq               = concatMap snd
                                 . map replaceCodon
                                 . zip (codonSplit germSeq)
                                 . codonSplit
    replaceCodon (x, y)
        | hamming x y > codonMut = ("---", "---")
        | otherwise              = (x, y)
    codonSplit                     = fullCodon . Split.chunksOf 3
    fullCodon                      = filter ((==3) . length)

-- Remove clone sequences that have stop codons in the first stopRange
-- codons
removeStopsCloneMap :: GeneticUnit -> Int -> CloneMap -> CloneMap
removeStopsCloneMap genUnit stopRange = M.map (filter (filterStops genUnit))
  where
    filterStops Nucleotide = not
                           . elem '*'
                           . take stopRange
                           . translate
                           . fastaSeq
    filterStops AminoAcid  = not . elem '*' . take stopRange . fastaSeq

-- Remove sequences that do not contain the string customFilter in the
-- customField location, split by "|". Note that this is 1 indexed and
-- 0 means to search the entire header for the customFilter.
removeCustomFilter :: Int -> String -> CloneMap -> CloneMap
removeCustomFilter 0 customFilter = M.map (filter inField)
  where
    inField = isInfixOf customFilter . fastaInfo
removeCustomFilter customField customFilter = M.map (filter inCustomField)
  where
    inCustomField x = isInfixOf customFilter
                    $ (Split.splitOn "|" . fastaInfo $ x) !! (customField - 1)

-- Remove clones that do not have any sequences after the filtrations
removeEmptyClone :: CloneMap -> CloneMap
removeEmptyClone = M.filter (not . null)

