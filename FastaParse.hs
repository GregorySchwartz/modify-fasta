-- FastaParse module.
-- By G.W. Schwartz
--
-- Collection of functions for the parsing of a fasta file.

module FastaParse where

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types

-- Parse a fasta file to return a list of FastaSequences
fastaParser :: String -> [FastaSequence]
fastaParser = map makeFastaSequence . Split.splitOn ">"
  where
    makeFastaSequence x = FastaSequence { fastaInfo = head . lines $ x
                                        , fastaSeq  = getSeq x
                                        }
    getSeq              = concat . drop 1 . lines

filterFasta :: [FastaSequence] -> [FastaSequence]
filterFasta = filter (\x -> not . elem '>' . fastaInfo $ x)
