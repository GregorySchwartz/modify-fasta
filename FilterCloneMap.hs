-- FilterCloneMap module.
-- By G.W. Schwartz
--
-- Collection of functions for the filtering of a CloneMap

module FilterCloneMap where

-- Built in
import Data.List
import Data.Char
import Data.Maybe
import Control.Monad.State
import qualified Data.Map as M

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types
import Diversity
import Translation

-- Remove highly mutated sequences (sequences with more than a third of
-- their sequence being mutated).
filterHighlyMutated :: GeneticUnit -> CloneMap -> CloneMap
filterHighlyMutated genUnit = M.mapWithKey filterMutated
  where
    filterMutated k xs  = filter (not . isHighlyMutated (snd k)) xs
    isHighlyMutated k x = ( (genericLength (readSeq genUnit . fastaSeq $ x)
                           :: Double) / 3 )
                       <= ( ( genericLength
                            . realMutations (readSeq genUnit . fastaSeq $ k)
                            $ readSeq genUnit . fastaSeq $ x) :: Double )
    realMutations k x   = filterCodonMutStab (\(x, y) -> x /= y)
                        . map snd
                        . mutation k
                        $ x
    filterCodonMutStab isWhat = filter (filterRules isWhat)
    filterRules isWhat x = isWhat x
                        && not (inTuple '-' x)
                        && not (inTuple '.' x)
                        && not (inTuple '~' x)
                        && not (inTuple 'N' x)
    inTuple c (x, y)
        | c == x || c == y = True
        | otherwise        = False
    mutation x y        = zip [1..] . zip x $ y
    readSeq Nucleotide  = id
    readSeq AminoAcid   = translate

-- Replace codons that have more than CodonMut mutations (make them "---"
-- codons).
removeCodonMutCount :: CodonMut -> String -> String -> CloneMap -> CloneMap
removeCodonMutCount codonMut codonMutType mutType = M.mapWithKey mapRemove
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
        | (codonMutOp codonMutType) (hamming x y) codonMut
       && isMutType (map toUpper mutType) x y          = (x, y)
        | otherwise                                    = ("---", "---")
    codonSplit                   = fullCodon . Split.chunksOf 3
    fullCodon                    = filter ((==3) . length)
    codonMutOp ">" = (>)
    codonMutOp "<" = (<)
    codonMutOp "=" = (==)
    isMutType "REPLACEMENT" x y = codon2aa x /= codon2aa y
    isMutType "SILENT" x y      = codon2aa x == codon2aa y
    isMutType _ x y             = True

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
-- 0 means to search the entire header for the customFilter. If the
-- customRemove option is enabled, this function will instead remove
-- sequences that have headers which match the custom filter, as opposed to
-- the other way around (this is defined in the "equal" function). Also
-- takes into account whether to filter on the germline versus the actual
-- sequences.
removeCustomFilter :: Bool
                   -> Bool
                   -> Bool
                   -> Maybe Int
                   -> String
                   -> State CloneMap ()
removeCustomFilter germ rm infixField customField customFilter
    | germ && ((customField == Just 0) || (isNothing customField))
        = state $ \cm -> ((), M.filterWithKey (\(p, k) _ -> inField k) cm)
    | germ && customField > Just 0
        = state $ \cm -> ((), M.filterWithKey (\(p, k) _ -> inCustomField k) cm)
    | (customField == Just 0) || (isNothing customField)
        = state $ \cm -> ((), M.map (filter inField) cm)
    | customField > Just 0 =
        state $ \cm -> ((), M.map (filter inCustomField) cm)
  where
    inField         = equal rm infixField customFilter . fastaInfo
    inCustomField x = equal rm infixField customFilter
                    . (!!) (Split.splitOn "|" . fastaInfo $ x)
                    $ (fromJust customField - 1)
    equal False False x = (==) x
    equal False True x  = isInfixOf x
    equal True False x  = (/=) x
    equal True True x   = not . isInfixOf x

removeAllCustomFilters :: Bool
                       -> Bool
                       -> Bool
                       -> [(Maybe Int, String)]
                       -> State CloneMap ()
removeAllCustomFilters germ rm infixField customFilters = do
    forM_ customFilters ( \(x, y) -> do
                          removeCustomFilter germ rm infixField x y )

-- -- Remove sequences that do not contain the string customFilter in the
-- -- customField location, split by "|". Note that this is 1 indexed and
-- -- 0 means to search the entire header for the customFilter. If the
-- -- customRemove option is enabled, this function will instead remove
-- -- sequences that have headers which match the custom filter, as opposed to
-- -- the other way around (this is defined in the "equal" function). Also
-- -- takes into account whether to filter on the germline versus the actual
-- -- sequences.
-- removeCustomFilter :: Bool
--                    -> Bool
--                    -> Bool
--                    -> Int
--                    -> String
--                    -> CloneMap
--                    -> State CloneMap ()
-- removeCustomFilter germ rm infixField customField customFilter
--     | germ && customField == 0 = state $ M.filterWithKey (\(p, k) _ -> inField k)
--     | germ && customField > 0  = M.filterWithKey (\(p, k) _ -> inCustomField k)
--     | customField == 0 = M.map (filter inField)
--     | customField /= 0 = M.map (filter inCustomField)
--   where
--     inField         = equal rm infixField customFilter . fastaInfo
--     inCustomField x = equal rm infixField customFilter
--                     $ (Split.splitOn "|" . fastaInfo $ x) !! (customField - 1)
--     equal False False x = (==) x
--     equal False True x  = isInfixOf x
--     equal True False x  = (/=) x
--     equal True True x   = not . isInfixOf x

-- Remove clones that do not have any sequences after the filtrations
removeEmptyClone :: CloneMap -> CloneMap
removeEmptyClone = M.filter (not . null)

-- Convert sequences to amino acids
convertToAminoAcidsCloneMap :: CloneMap -> CloneMap
convertToAminoAcidsCloneMap = M.mapKeys keyMap
                            . M.map (map translateFastaSequence)
  where
    keyMap (x, y) = (x, translateFastaSequence y)
    translateFastaSequence x = x { fastaSeq = translate . fastaSeq $ x }
