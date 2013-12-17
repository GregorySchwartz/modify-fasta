-- Print module
-- By G.W. Schwartz
--
-- Collection of functions for the printing of data (converting data
-- structures into strings for use with writing to output files).

module Print where

-- Local
import Types

-- Return the results of the diversity analysis in string form for saving
-- to a file
printFasta :: [FastaSequence] -> String
printFasta fastaList = body
  where
    body             = unlines . map mapLine $ fastaList
    mapLine x = ">" ++ fastaInfo x ++ "\n" ++ fastaSeq x
