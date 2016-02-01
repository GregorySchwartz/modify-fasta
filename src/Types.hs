-- Types module.
-- By Gregory W. Schwartz
--
-- Collects all application specific types.

module Types where

-- Built-in
import qualified Data.Text as T
import qualified Data.Map.Strict as Map

-- Cabal
import Data.Fasta.Text.Types

-- Algebraic
data GeneticUnit   = AminoAcid | Nucleotide deriving (Read, Show)
data FrameType     = InFrame | OutFrame deriving (Read, Show)

-- Basic
type ID       = Int
type Codon    = T.Text
type CodonMut = Int
type Field    = Int
type Start    = Int
type Stop     = Int
type Position = Int
type Frame    = Int

-- Advanced
type CloneEntry     = (Germline, [FastaSequence])
type FillInValue    = (Field, Start, Char)
type CodonTable     = [(T.Text, Char)]
type Mutation       = (Char, Char)
type CountMap       = Map.Map (Position, Mutation) Int
type CodonMutations = [[(Position, Mutation)]]
