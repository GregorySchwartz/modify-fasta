-- FilterCloneMap module.
-- By G.W. Schwartz
--
-- Collection of functions for the filtering of a CloneMap

{-# LANGUAGE BangPatterns #-}

module FilterCloneMap where

-- Built in
import Data.List
import Data.Char
import Data.Maybe
import Data.Either
import qualified Data.Set as S
import qualified Data.Map as M
import Text.Regex.TDFA

-- Cabal
import qualified Data.List.Split as Split
import Data.Fasta.String

-- Local
import Types
import Diversity

-- Altered version of listToMaybe
listToMaybe' :: [a] -> Maybe [a]
listToMaybe' []      = Nothing
listToMaybe' x       = Just x

-- Returns the string of messed up sequences after trying to translate
getErrorString :: CloneMap -> Maybe String
getErrorString = listToMaybe'
               . unlines
               . filter (not . null)
               . map snd
               . M.toAscList
               . M.map (intercalate "\n" . lefts . map translate)

-- Remove highly mutated sequences (sequences with more than a third of
-- their sequence being mutated).
filterHighlyMutated :: GeneticUnit -> CloneMap -> (CloneMap, (Maybe String))
filterHighlyMutated !genUnit !cloneMap = (newCloneMap, errorString)
  where
    newCloneMap         = M.map (map snd . filter fst . rights) errorCloneMap
    errorString         = listToMaybe'
                        . unlines
                        . filter (not . null)
                        . map snd
                        . M.toAscList
                        . M.map (intercalate "\n" . lefts)
                        $ errorCloneMap
    errorCloneMap       = M.mapWithKey assignMutated cloneMap
    assignMutated k xs  = map (isHighlyMutated (snd k)) xs
    isHighlyMutated k x =
        case (readSeq genUnit k, readSeq genUnit x) of
            ((Right a), (Right b)) -> (\n -> Right (n, b))
                                    $ ( (genericLength (fastaSeq a) :: Double)
                                      / 3 )
                                   <= ( ( genericLength
                                        . realMutations (fastaSeq a)
                                        $ fastaSeq b )
                                       :: Double )
            ((Left a), (Right _)) -> Left (unwords ["Germline: ", a])
            ((Right _), (Left b)) -> Left (unwords ["Sequence: ", b])
            ((Left a), (Left b))  -> Left ( unwords [ "Sequence:"
                                                    , b
                                                    , "with Germline:"
                                                    , a ] )
    realMutations k x   = filterCodonMutStab (\(y, z) -> y /= z)
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
    mutation x y        = zip [1..] . zip x $ y
    readSeq Nucleotide x = Right (id x)
    readSeq AminoAcid x  = translate x

-- Replace codons that have more than CodonMut mutations (make them "---"
-- codons).
removeCodonMutCount :: CodonMut -> String -> String -> CloneMap -> CloneMap
removeCodonMutCount codonMut codonMutType mutType = M.mapWithKey mapRemove
  where
    mapRemove (_, germ) clones   = map (removeCodon germ) clones
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
    isMutType _ _ _             = True

-- Remove clone sequences that have stop codons in the first stopRange
-- codons
removeStopsCloneMap :: GeneticUnit
                    -> Int
                    -> CloneMap
                    -> (CloneMap, Maybe String)
removeStopsCloneMap !genUnit !stopRange !cloneMap = ( newCloneMap
                                                    , errorString )
  where
    errorString = getErrorString cloneMap
    newCloneMap = M.map (filter (filterStops genUnit)) cloneMap
    filterStops Nucleotide x = (isRight' . translate $ x)
                            && ( not
                               . elem '*'
                               . take stopRange
                               . fastaSeq
                               . fromEither
                               . translate ) x
    filterStops AminoAcid  x = not . elem '*' . take stopRange . fastaSeq $ x
    isRight' (Right _)       = True
    isRight' _               = False
    fromEither (Right x)     = x
    fromEither (Left _)      = error "This should not have happened"

-- Remove duplicate sequences
removeDuplicatesCloneMap :: CloneMap -> CloneMap
removeDuplicatesCloneMap cloneMap = M.map
                                    (filter (\x -> S.member x duplicateSet))
                                  $ cloneMap
  where
    duplicateSet = S.fromList
                 . nubBy (\x y -> fastaSeq x == fastaSeq y)
                 . concatMap snd
                 . M.toAscList
                 $ cloneMap

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
                   -> Maybe Int
                   -> String
                   -> CloneMap
                   -> CloneMap
removeCustomFilter germ rm customField customFilter cloneMap
    | germ && ((customField == Just 0) || (isNothing customField))
        = M.filterWithKey (\(_, k) _ -> inField k) cloneMap
    | germ && customField > Just 0
        = M.filterWithKey (\(_, k) _ -> inCustomField k) cloneMap
    | (customField == Just 0) || (isNothing customField)
        = M.map (filter inField) cloneMap
    | customField > Just 0 =
        M.map (filter inCustomField) cloneMap
  where
    inField         = equal rm customFilter . fastaHeader
    inCustomField x = equal rm customFilter
                    . (!!) (Split.splitOn "|" . fastaHeader $ x)
                    $ (fromJust customField - 1)
    equal False x y = y =~ x :: Bool
    equal True x y  = not . equal False x $ y

removeAllCustomFilters :: Bool
                       -> Bool
                       -> CloneMap
                       -> [(Maybe Int, String)]
                       -> CloneMap
removeAllCustomFilters germ rm cloneMap = foldl' filterMap cloneMap
  where
    filterMap acc (x, y) = removeCustomFilter germ rm x y acc

-- Remove clones that do not have any sequences after the filtrations
removeEmptyClone :: CloneMap -> CloneMap
removeEmptyClone = M.filter (not . null)

-- Convert sequences to amino acids
convertToAminoAcidsCloneMap :: CloneMap -> (CloneMap, (Maybe String))
convertToAminoAcidsCloneMap !cloneMap = (newCloneMap, errorString)
  where
    errorString = getErrorString cloneMap
    newCloneMap = M.mapKeys keyMap
                . M.map (map translateFastaSequence)
                $ cloneMap
    keyMap (x, y) = (x, translateFastaSequence y)
    translateFastaSequence x = x { fastaSeq = fastaSeq
                                            . fromEither
                                            . translate
                                            $ x }
    fromEither (Right x)     = x
    fromEither (Left _)      = error "This should not have happened"
