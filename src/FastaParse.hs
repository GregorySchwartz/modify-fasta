-- FastaParse module.
-- By G.W. Schwartz
--
-- Collection of functions for the parsing of a fasta file.

module FastaParse where

-- Built in
import Data.List
import Data.Char
import qualified Data.Map as M

-- Cabal
import qualified Data.List.Split as Split
import Data.Fasta.String

-- Takes a clip fasta file string and removes newlines in the sequences to
-- make this compatible with the fasta parser. The lineCompress function
-- should get rid of any extra newlines messing with the code. Also a flag
-- to see if we should remove N's or not.
joinSeq :: Bool -> String -> String
joinSeq removeNFlag = lineCompress
                    . tail
                    . concat
                    . map newEntry
                    . filter (/= "")
                    . Split.splitOn ">>"
  where
    newEntry x             = if elem '>' x then cloneEntry x else germEntry x
    germEntry x            = newGerm x
    cloneEntry x           = newGerm (germline x)
                          ++ concat (map newClone . filter (/= "") . clone $ x)
    newGerm x
        | getSeq x /= ""   = "\n>>" ++ (header x) ++ "\n" ++ (getSeq x)
        | otherwise        = ""
    newClone x
        | getSeq x /= ""   = "\n>" ++ (header x) ++ "\n" ++ (getSeq x)
        | otherwise        = ""
    germline               = head . Split.splitOn ">"
    clone                  = tail . Split.splitOn ">"
    header                 = head . lines
    getSeq                 = map (removeN removeNFlag) . concat . tail . lines
    lineCompress []        = []
    lineCompress ('\n':xs) = '\n' : (lineCompress $ dropWhile (== '\n') xs)
    lineCompress (x:xs)    = x : (lineCompress xs)
    removeN False x        = x
    removeN True 'N'       = '-'
    removeN True 'n'       = '-'
    removeN True x         = x

-- Takes a fasta file string of the format ">>[Germline header]\n[Germline
-- sequence]\n>[Mutant header]\n[Mutant
-- sequence]\n>[MH]\n[MS]...\n>>[GH]\n[GS]\n>[MH]\n[MS]...etc" and returns
-- a CloneMap in order to generate the basic building block for the
-- mutation counting.  Note: Several repeating germlines, so they need
-- a unique identifier (an integer in this case). Also removes empty clones
-- (only a germline appears).
generateCloneMap :: String -> CloneMap
generateCloneMap = M.fromList . assocList
  where
    assocList                      = map assocMap . germlineSplit
    assocMap (x, y)                = ((x, germline y), clones y)
    germlineSplit                  = zip [0..]
                                   . filter (\x -> elem '>' x)
                                   . filter (/= "")
                                   . Split.splitOn ">>"
    germline x                     = FastaSequence
                                   { fastaHeader = head . lines $ x
                                   , fastaSeq  = map toUpper
                                               . head
                                               . drop 1
                                               . lines
                                               $ x
                                   }
    clones                         = clones2FastaSequence . drop 2 . lines
    clones2FastaSequence []            = []
    clones2FastaSequence (h:s:xs) = FastaSequence { fastaHeader = h
                                                  , fastaSeq    = s
                                                  }
                                  : (clones2FastaSequence xs)

-- Adds filler germlines to normal fasta files
addFillerGermlines :: String -> String
addFillerGermlines = replace ">" ">>filler\n---\n>"
  where
    replace old new = intercalate new . Split.splitOn old
