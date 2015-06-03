-- Utility module
-- By Gregory W. Schwartz
--
-- Collects utility functions for the main files

module Utility ( addLengthHeader
               , addFillerGermlines ) where

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
