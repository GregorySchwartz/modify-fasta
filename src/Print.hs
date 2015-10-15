-- Print module
-- By Gregory W. Schwartz
--
-- Collection of functions for the printing of data (converting data
-- structures into strings for use with writing to output files).

{-# LANGUAGE BangPatterns, OverloadedStrings #-}

module Print where

-- Built-in
import Data.List
import qualified Data.Map as M
import qualified Data.Text as T

-- Cabal
import Data.Fasta.Text
import TextShow

-- Local
import Types

-- Return the results of the filtration in text form for saving
-- to a file
printFasta :: CloneMap -> T.Text
printFasta = body
  where
    body                = T.unlines
                        . map mapGerm
                        . M.toAscList
                        . M.map (T.intercalate "\n" . map showFasta)
    mapGerm ((_, y), z) = mconcat [ ">>"
                                  , fastaHeader y
                                  , "\n"
                                  , fastaSeq y
                                  , "\n"
                                  , z
                                  ]

-- Return the results of the filtration in text form for saving
-- to a file and excluding germline
printFastaNoGermline :: CloneMap -> T.Text
printFastaNoGermline = body
  where
    body                = T.unlines
                        . map mapGerm
                        . M.toAscList
                        . M.map (T.intercalate "\n" . map showFasta)
    mapGerm ((_, _), z) = z

printSequenceCount :: Bool -> Int -> CloneMap -> T.Text
printSequenceCount clip idx s = body
  where
    body = T.unlines [ "Allele List: "
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
                     , mappend "Number of genes: " (showt . M.size $ germlineMap)
                     , mappend "Number of alleles: " (showt . M.size $ alleleMap)
                     ]
    geneCounts     = T.intercalate "\n" . map toLine . M.toAscList $ germlineMap
    alleleCounts   = T.intercalate "\n" . map fst . M.toAscList $ alleleMap
    toLine (x, y)  = x
           `mappend` (T.replicate (30 - T.length x) " ")
           `mappend` showt y
    germlineMap    = M.fromListWith (+)
                   . map (\(x, y) -> (head . T.splitOn "*" $ x, y))
                   $ geneAlleleList
    alleleMap      = M.fromListWith (+) geneAlleleList
    geneAlleleList = map (countProp clip) . M.toAscList $ s
    countProp True ((_, x), y)  = (getField idx x, length y)
    countProp False ((_, _), y) = (getField idx . head $ y, 1)
    getField f h   = splitHeader h !! (f - 1)
    splitHeader    = T.splitOn "|" . fastaHeader

-- | Takes a clone entry and returns a formatted text with or without
-- germline
printCloneEntry :: Bool -> CloneEntry -> T.Text
printCloneEntry False (!germline, !fseqs)  =
    T.intercalate "\n" [ T.cons '>' $ showFasta germline
                       , T.intercalate "\n" . map showFasta $ fseqs
                       ]
printCloneEntry True (!germline, !fseqs) = T.intercalate "\n"
                                          . map showFasta
                                          $ fseqs
