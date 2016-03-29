-- TransformFastaList module.
-- By Gregory W. Schwartz
--
-- Collection of functions that transform a fasta sequence in some way

{-# LANGUAGE BangPatterns, OverloadedStrings #-}

module TransformFastaList ( convertToAminoAcidsFastaSequence
                          , replaceChars
                          , fillInSequence
                          , changeField
                          , changeAllFields
                          , getRegionSequence
                          , trimFasta
                          , removeUnknownNucs
                          ) where

-- Built-in
import Data.List
import Data.Char
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text
import Text.Regex.TDFA
import Text.Regex.TDFA.Text

-- Local
import Types
import Utility

-- | Convert sequences to amino acids
convertToAminoAcidsFastaSequence :: CodonTable
                                 -> FastaSequence
                                 -> FastaSequence
convertToAminoAcidsFastaSequence table = fromEither . customTranslate table 1

-- | Fill in the sequence with corrected nucleotides or amino acids
fillInSequence :: Field -> Start -> Char -> FastaSequence -> FastaSequence
fillInSequence f s c fs = fs { fastaSeq = newFastaSeq }
  where
    newFastaSeq  = (\xs -> if T.singleton c `T.isInfixOf` xs then "" else xs)
                 $ first `mappend` replaceChars c old new
    new          = (T.splitOn "|" . fastaHeader $ fs) !! (f - 1)
    (first, old) = T.splitAt (s - 1) . fastaSeq $ fs

-- | Change a field to a match, so a regex "ch.*_" to field 2 of
-- ">abc|brie_cheese_dude" would result in ">abc|cheese_". Useful for
-- getting specific properties from a field
changeField :: Maybe Field -> T.Text -> FastaSequence -> FastaSequence
changeField Nothing _ fs          = fs
changeField (Just field) regex fs = fs { fastaHeader = newFastaHeader }
  where
    newFastaHeader  = T.intercalate "|"
                    . F.toList
                    . Seq.update (field - 1) newField
                    $ splitField
    newField        = (Seq.index splitField (field - 1)) =~ regex :: T.Text
    splitField      = Seq.fromList . T.splitOn "|" . fastaHeader $ fs

-- | Change all fields to their matches based on changeField
changeAllFields :: FastaSequence -> [(Maybe Int, T.Text)] -> FastaSequence
changeAllFields = foldl' (\fs (!x, !y) -> changeField x y fs)

-- | Get a region of a text, 0 indexed
getRegion :: Maybe Start -> Maybe Stop -> T.Text -> T.Text
getRegion Nothing Nothing          = id
getRegion Nothing (Just stop)      = T.take stop
getRegion (Just start) Nothing     = T.drop start
getRegion (Just start) (Just stop) = T.take (stop - start) . T.drop start

-- | Get a region of a sequence, 1 indexed
getRegionSequence :: Maybe Start -> Maybe Stop -> FastaSequence -> FastaSequence
getRegionSequence start0 stop fs = fs { fastaSeq = newFastaSeq }
  where
    newFastaSeq = getRegion start stop . fastaSeq $ fs
    start       = fmap (flip (-) 1) start0

-- | Trim the sequence. For default, trim the extra nucleotides off the
-- end. Otherwise, use UCSC frames (0, 1, 2) to cut off the beginning (0
-- means in frame) or the end (0 means the nucleotide AFTER the end of the
-- sequence is in frame). For amino acids, if it's not in frame, just cut
-- that amino acid.
trim :: GeneticUnit -> Maybe FrameType -> Maybe Frame -> T.Text -> T.Text
trim Nucleotide Nothing Nothing x              = T.dropEnd
                                                 (T.length x `mod` 3)
                                                 x
trim _ _ Nothing x                             = x
trim _ _ (Just 0) x                            = x
trim Nucleotide (Just InFrame) (Just frame) x  = T.drop (3 - frame) x
trim Nucleotide (Just OutFrame) (Just frame) x = T.dropEnd frame x
trim AminoAcid (Just InFrame) (Just _) x       = T.drop 1 x
trim AminoAcid (Just OutFrame) (Just _) x      = T.dropEnd 1 x

-- | Trim off extra nucleotides (or amino acids) from a fasta sequence. If
-- inframe and outframe are specified, instead cut off based on those frames.
trimFasta :: GeneticUnit
          -> Maybe Frame
          -> Maybe Frame
          -> FastaSequence
          -> FastaSequence
trimFasta Nucleotide Nothing Nothing fs =
    fs { fastaSeq = trim Nucleotide Nothing Nothing . fastaSeq $ fs }
trimFasta genU inF outF fs = fs { fastaSeq = trim genU (Just InFrame) inF
                                           . trim genU (Just OutFrame) outF
                                           . fastaSeq $ fs
                                }

-- | Convert non standard nucleotides to gaps
removeUnknownNucs :: FastaSequence -> FastaSequence
removeUnknownNucs fs = fs { fastaSeq = T.map changeNuc . fastaSeq $ fs }
  where
    changeNuc x
        | toUpper x `elem` ("ATCGN.-" :: String) = x
        | otherwise                              = '-'

