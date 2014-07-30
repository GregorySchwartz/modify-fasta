-- Print module
-- By G.W. Schwartz
--
-- Collection of functions for the printing of data (converting data
-- structures into strings for use with writing to output files).

module Print where

-- Built-in
import Data.List
import qualified Data.Map as M

-- Cabal
import qualified Data.List.Split as Split

-- Local
import Types

-- Return the results of the filtration in string form for saving
-- to a file
printFasta :: CloneMap -> String
printFasta = body
  where
    body                = unlines
                        . map mapGerm
                        . M.toAscList
                        . M.map (intercalate "\n" . map mapClone)
    mapGerm ((_, y), z) = ">>" ++ fastaInfo y ++ "\n" ++ fastaSeq y ++ "\n" ++ z
    mapClone x          = fastaInfo x ++ "\n" ++ fastaSeq x

-- Return the results of the filtration in string form for saving
-- to a file and excluding germline
printFastaNoGermline :: CloneMap -> String
printFastaNoGermline = body
  where
    body                = unlines
                        . map mapGerm
                        . M.toAscList
                        . M.map (intercalate "\n" . map mapClone)
    mapGerm ((_, _), z) = z
    mapClone x          = fastaInfo x ++ "\n" ++ fastaSeq x

printSequenceCount :: Bool -> Int -> CloneMap -> String
printSequenceCount clip idx s = body
  where
    body = unlines [ "Allele List: "
                   , ""
                   , "----------------------------------------------------"
                   , ""
                   , alleleCounts
                   , ""
                   , "----------------------------------------------------"
                   , ""
                   , "Gene List: "
                   , ""
                   , "----------------------------------------------------"
                   , ""
                   , geneCounts
                   , ""
                   , "----------------------------------------------------"
                   , ""
                   , "Number of genes: " ++ (show . M.size $ germlineMap)
                   , "Number of alleles: " ++ (show . M.size $ alleleMap)
                   ]
    geneCounts     = intercalate "\n" . map toLine . M.toAscList $ germlineMap
    alleleCounts   = intercalate "\n" . map fst . M.toAscList $ alleleMap
    toLine (x, y)  = x ++ (take (30 - length x) . repeat $ ' ') ++ show y
    germlineMap    = M.fromListWith (+)
                   . map (\(x, y) -> (head . Split.splitOn "*" $ x, y))
                   $ geneAlleleList
    alleleMap      = M.fromListWith (+) geneAlleleList
    geneAlleleList = map (countProp clip) . M.toAscList $ s
    countProp True ((_, x), y)  = (getField idx x, length y)
    countProp False ((_, _), y) = (getField idx . head $ y, 1)
    getField f h   = splitHeader h !! (f - 1)
    splitHeader    = Split.splitOn "|" . fastaInfo

