-- Types module.
-- By Gregory W. Schwartz
--
-- Collects all application specific types.

module Types where

-- Built-in
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text.Types

-- Algebraic
data GeneticUnit   = AminoAcid | Nucleotide deriving (Read, Show)

-- Basic
type ID       = Int
type Codon    = T.Text
type CodonMut = Int
type Field    = Int
type Start    = Int
type Stop     = Int

-- Advanced
type CloneEntry  = (Germline, [FastaSequence])
type FillInValue = (Field, Start, Char)
