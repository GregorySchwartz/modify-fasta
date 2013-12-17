-- Types module.
-- By G.W. Schwartz
--
-- Collects all application specific types.

module Types where

-- Algebraic
data FastaSequence = FastaSequence { fastaInfo :: String
                                   , fastaSeq  :: String
                                   } deriving (Eq, Ord, Show)
