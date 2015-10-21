-- FilterFastaList module.
-- By Gregory W. Schwartz
--
-- Collection of functions for the filtering of a pipesFasta

{-# LANGUAGE OverloadedStrings, FlexibleContexts #-}

module FilterFastaList ( hasNoStops
                       , isInFrame
                       , hasCustomFilter
                       , hasAllCustomFilters
                       ) where

-- Built in
import Data.List
import Data.Maybe
import Text.Regex.TDFA
import Text.Regex.TDFA.Text
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text

-- Local
import Types

-- Remove clone sequences that have stop codons in the first stopRange
-- codons
hasNoStops :: GeneticUnit
           -> Int
           -> FastaSequence
           -> Bool
hasNoStops genUnit stopRange = result . stop genUnit
  where
    result (Right x) = x
    result (Left x)  = error . T.unpack $ x
    stop Nucleotide = fmap ( not
                           . T.isInfixOf "*"
                           . T.take stopRange
                           . fastaSeq )
                    . translate 1
    stop AminoAcid = Right . not . T.isInfixOf "*" . T.take stopRange . fastaSeq

-- Remove out of frame sequences
isInFrame :: FastaSequence -> Bool
isInFrame = (== 0)
          . mod 3
          . T.length
          . T.filter (\x -> not . T.isInfixOf (T.singleton x) $ ".-")
          . fastaSeq

-- Remove sequences that do not contain the string customFilter in the
-- customField location, split by "|". Note that this is 1 indexed and
-- 0 means to search the entire header for the customFilter. If the
-- customRemove option is enabled, this function will instead remove
-- sequences that have headers which match the custom filter, as opposed to
-- the other way around (this is defined in the "equal" function). Also
-- takes into account whether to filter on the germline versus the actual
-- sequences.
hasCustomFilter :: Bool
                -> Maybe Int
                -> T.Text
                -> FastaSequence
                -> Bool
hasCustomFilter rm customField customFilter fasta
    | customField == Just 0 || isNothing customField = inField fasta
    | customField > Just 0                           = inCustomField fasta
  where
    inField         = equal rm customFilter . fastaHeader
    inCustomField x = equal rm customFilter
                    . (!!) (T.splitOn "|" . fastaHeader $ x)
                    $ (fromJust customField - 1)
    equal :: Bool -> T.Text -> T.Text -> Bool
    equal False x y = y =~ x :: Bool
    equal True x y  = not . equal False x $ y

hasAllCustomFilters :: Bool
                    -> [(Maybe Int, T.Text)]
                    -> FastaSequence
                    -> Bool
hasAllCustomFilters rm filters f = all filterMap filters
  where
    filterMap (x, y) = hasCustomFilter rm x y f
