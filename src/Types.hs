-- Types module.
-- By Gregory W. Schwartz
--
-- Collects all application specific types.

module Types where

-- Algebraic
data GeneticUnit   = AminoAcid | Nucleotide deriving (Read, Show)

-- Basic
type ID         = Int
type Codon      = String
type CodonMut   = Int
