-- FilterCloneMap module.
-- By G.W. Schwartz
--
-- Collection of functions for the filtering of a CloneMap

{-# LANGUAGE BangPatterns, OverloadedStrings, FlexibleContexts #-}

module FilterCloneMap where

-- Built in
import Data.List
import Data.Char
import Data.Maybe
import Data.Either
import qualified Data.Set as S
import qualified Data.Map as M
import Text.Regex.TDFA
import Text.Regex.TDFA.Text
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text

-- Local
import Types
import Diversity

-- | Check if the data structure is Right
isRight' :: Either a b -> Bool
isRight' (Right _)       = True
isRight' _               = False

-- | Altered version of listToMaybe
listToMaybe' :: [a] -> Maybe [a]
listToMaybe' []      = Nothing
listToMaybe' x       = Just x

-- | Remove highly mutated sequences (sequences with more than a third of
-- their sequence being mutated).
filterHighlyMutated :: GeneticUnit
                    -> CodonTable
                    -> CloneMap
                    -> (CloneMap, Maybe String)
filterHighlyMutated !genUnit !table !cloneMap = (newCloneMap, errorString)
  where
    newCloneMap           = M.map (map snd . filter (not . fst) . rights)
                            errorCloneMap
    errorString           = listToMaybe'
                          . unlines
                          . filter (not . null)
                          . map snd
                          . M.toAscList
                          . M.map (intercalate "\n" . lefts)
                          $ errorCloneMap
    errorCloneMap         = M.mapWithKey assignMutated cloneMap
    assignMutated k       = map (isHighlyMutated (snd k))
    isHighlyMutated !k !x =
        case (readSeq genUnit k, readSeq genUnit x) of
            ((Right a), (Right b)) -> (\n -> Right (n, b))
                                    $ ( (fromIntegral (T.length (fastaSeq a)) :: Double)
                                      / 3 )
                                   <= ( ( genericLength
                                        . realMutations (fastaSeq a)
                                        $ fastaSeq b ) )
            ((Left a), (Right _)) -> Left (unwords ["Germline: ", T.unpack a])
            ((Right _), (Left b)) -> Left (unwords ["Sequence: ", T.unpack b])
            ((Left a), (Left b))  -> Left (unwords [ "Sequence:"
                                                    , T.unpack b
                                                    , "with Germline:"
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

-- | Replace codons that have more than CodonMut mutations (make them "---"
-- codons).
removeCodonMutCount :: CodonMut -> T.Text -> T.Text -> CloneMap -> CloneMap
removeCodonMutCount codonMut codonMutType mutType = M.mapWithKey mapRemove
  where
    mapRemove (_, germ)          = map (removeCodon germ)
    removeCodon germ clone       = clone { fastaSeq
                                         = remove (fastaSeq germ)
                                         . fastaSeq $ clone }
    remove germSeq               = mconcat
                                 . map (snd . replaceCodon)
                                 . zip (codonSplit germSeq)
                                 . codonSplit
    replaceCodon (x, y)
        | (codonMutOp codonMutType) (hamming x y) codonMut
       && isMutType (T.toUpper mutType) x y          = (x, y)
        | otherwise                                    = ("---", "---")
    codonSplit                   = fullCodon . T.chunksOf 3
    fullCodon                    = filter ((== 3) . T.length)
    codonMutOp ">" = (>)
    codonMutOp "<" = (<)
    codonMutOp "=" = (==)
    isMutType "REPLACEMENT" x y = codon2aa x /= codon2aa y
    isMutType "SILENT" x y      = codon2aa x == codon2aa y
    isMutType _ _ _             = True

-- | Remove clone sequences that have stop codons in the first stopRange
-- codons
removeStopsCloneMap :: GeneticUnit
                    -> CodonTable
                    -> Int
                    -> CloneMap
                    -> (CloneMap, Maybe String)
removeStopsCloneMap !genUnit !table !stopRange !cloneMap = ( newCloneMap
                                                           , errorString )
  where
    errorString = listToMaybe'
                . unlines
                . filter (not . null)
                . map snd
                . M.toAscList
                . M.map ( intercalate "\n"
                        . map T.unpack
                        . lefts
                        . map (customTranslate table 1)
                        )
                $ cloneMap
    newCloneMap = M.map (filter (filterStops genUnit)) cloneMap
    filterStops Nucleotide x = (isRight' . customTranslate table 1 $ x)
                            && ( not
                               . T.isInfixOf "*"
                               . T.take stopRange
                               . fastaSeq
                               . fromEither
                               . customTranslate table 1 ) x
    filterStops AminoAcid  x = not
                             . T.isInfixOf "*"
                             . T.take stopRange
                             . fastaSeq
                             $ x
    fromEither (Right x)     = x
    fromEither (Left x)      = error (T.unpack x)

-- | Remove duplicate sequences
removeDuplicatesCloneMap :: CloneMap -> CloneMap
removeDuplicatesCloneMap cloneMap = M.map
                                    (filter (`S.member` duplicateSet))
                                    cloneMap
  where
    duplicateSet = S.fromList
                 . nubBy (\x y -> fastaSeq x == fastaSeq y)
                 . concatMap snd
                 . M.toAscList
                 $ cloneMap

-- | Remove out of frame sequences
removeOutOfFrameSeqs :: CloneMap -> CloneMap
removeOutOfFrameSeqs = M.map (filter isInFrame)
  where
    isInFrame  = (== 0)
               . mod 3
               . T.length
               . T.filter (\x -> not $ T.isInfixOf (T.singleton x) ".-")
               . fastaSeq

-- | Remove sequences that do not contain the string customFilter in the
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
                   -> T.Text
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
                    . (!!) (T.splitOn "|" . fastaHeader $ x)
                    $ (fromJust customField - 1)
    equal False x y = y =~ x :: Bool
    equal True x y  = not . equal False x $ y

removeAllCustomFilters :: Bool
                       -> Bool
                       -> CloneMap
                       -> [(Maybe Int, T.Text)]
                       -> CloneMap
removeAllCustomFilters germ rm = foldl' filterMap
  where
    filterMap acc (x, y) = removeCustomFilter germ rm x y acc

-- | Remove clones that do not have any sequences after the filtrations
removeEmptyClone :: CloneMap -> CloneMap
removeEmptyClone = M.filter (not . null)

-- | Convert sequences to amino acids
convertToAminoAcidsCloneMap :: CodonTable
                            -> CloneMap
                            -> (CloneMap, Maybe String)
convertToAminoAcidsCloneMap table cloneMap = (newCloneMap, errorString)
  where
    newCloneMap   = M.mapKeysWith (++) (\(!x, !y) -> (x, fromEither y))
                  . M.filterWithKey (\(_, !y) _ -> isRight' y)
                  . M.map rights
                  $ errorCloneMap
    errorString   = listToMaybe'
                  . concatMap snd
                  . M.toAscList
                  . M.mapWithKey (\(_, !y) v -> (++) (eitherToString y)
                                             . concatMap T.unpack
                                             . lefts
                                             $ v )
                  $ errorCloneMap
    errorCloneMap = M.mapKeys keyMap
                  . M.map (map (customTranslate table 1))
                  $ cloneMap
    keyMap (!x, !y) = (x, customTranslate table 1 y)
    eitherToString (Right _) = ""
    eitherToString (Left x)  = T.unpack x
    fromEither (Right x)     = x
    fromEither (Left x)      = error (T.unpack x)
