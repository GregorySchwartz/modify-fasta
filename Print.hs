-- Print module
-- By G.W. Schwartz
--
-- Collection of functions for the printing of data (converting data
-- structures into strings for use with writing to output files).

module Print where

-- Built-in
import Data.List
import qualified Data.Map as M

-- Local
import Types

-- Return the results of the filtration in string form for saving
-- to a file
printFasta :: CloneMap -> String
printFasta = body
  where
    body                = intercalate "\n"
                        . map mapGerm
                        . M.toAscList
                        . M.map (intercalate "\n" . map mapClone)
    mapGerm ((x, y), z) = ">>" ++ fastaInfo y ++ "\n" ++ fastaSeq y ++ "\n" ++ z
    mapClone x          = fastaInfo x ++ "\n" ++ fastaSeq x

-- Return the results of the filtration in string form for saving
-- to a file and excluding germline
printFastaNoGermline :: CloneMap -> String
printFastaNoGermline = body
  where
    body                = intercalate "\n"
                        . map mapGerm
                        . M.toAscList
                        . M.map (intercalate "\n" . map mapClone)
    mapGerm ((x, y), z) = z
    mapClone x          = fastaInfo x ++ "\n" ++ fastaSeq x
