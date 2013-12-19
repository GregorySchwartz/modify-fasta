-- Filter CLIP Fasta
-- By G.W. Schwartz

-- Takes a CLIP fasta file and filters the fasta file in several optional
-- ways.

-- Cabal
import Options.Applicative

-- Local
import Types
import FastaParse
import FilterCloneMap
import Print

-- Command line arguments
data Options = Options { input           :: String
                       , removeGermlines :: Bool
                       , removeStops     :: Bool
                       , inputCodonMut   :: CodonMut
                       , output          :: String
                       }

-- Command line options
options :: Parser Options
options = Options
      <$> strOption
          ( long "input"
         <> short 'i'
         <> metavar "FILE"
         <> value ""
         <> help "The input CLIP fasta file" )
      <*> switch
          ( long "removeGermlines"
         <> short 'g'
         <> help "Whether to remove germlines" )
      <*> switch
          ( long "removeStops"
         <> short 's'
         <> help "Whether to remove clone sequences with stop codons" )
      <*> option
          ( long "inputCodonMut"
         <> short 'c'
         <> metavar "[0]|1|2|3"
         <> value 0
         <> help "Only include codons with this many mutations or more\
                 \ (0 is the same as include all codons). Converts\
                 \ the codon to " )
      <*> strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> value "output.fasta"
         <> help "The output fasta file" )

filterCLIPFasta :: Options -> IO ()
filterCLIPFasta opts = do
    unfilteredContents <- readFile . input $ opts
    let contents              = joinSeq unfilteredContents

    let codonMut              = inputCodonMut opts

    let cloneMap              = generateCloneMap contents
    let cloneMapNoStops       = if (removeStops opts)
                                    then removeStopsCloneMap cloneMap
                                    else cloneMap
    let cloneMapNoEmptyClones = removeEmptyClone cloneMapNoStops

    let outputString = if (removeGermlines opts)
                           then  printFastaNoGermline cloneMapNoStops
                           else  printFasta cloneMapNoStops

    writeFile (output opts) outputString

main :: IO ()
main = execParser opts >>= filterCLIPFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Filter CLIP fasta files in several optional ways"
     <> header "Filter CLIP Fasta, Gregory W. Schwartz" )
