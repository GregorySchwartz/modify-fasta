-- Types module.
-- By G.W. Schwartz
--
-- Collects all application specific types.

module Types where

-- Built in
import qualified Data.Map as M

-- Algebraic
data FastaSequence = FastaSequence { fastaInfo :: String
                                   , fastaSeq  :: String
                                   } deriving (Eq, Ord, Show)

-- Basic
type ID         = Int
type Codon      = String
type Clone      = FastaSequence
type Germline   = FastaSequence
type CodonMut   = Int

-- Advanced
type CloneMap      = M.Map (ID, Germline) [Clone]
