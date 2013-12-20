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
                       , aminoAcids      :: GeneticUnit
                       , removeN         :: Bool
                       , removeGermlines :: Bool
                       , removeStops     :: Bool
                       , inputStopRange  :: Int
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
      <*> flag Nucleotide AminoAcid
          ( long "aminoAcids"
         <> short 'a'
         <> help "Whether these sequences are composed of\
                 \ amino acids (AminoAcid) or nucleotides (Nucleotide)" )
      <*> switch
          ( long "removeN"
         <> short 'N'
         <> help "Whether to remove N or n in the sequence" )
      <*> switch
          ( long "removeGermlines"
         <> short 'g'
         <> help "Whether to remove germlines" )
      <*> switch
          ( long "removeStops"
         <> short 's'
         <> help "Whether to remove clone sequences with stop codons" )
      <*> option
          ( long "inputStopRange"
         <> short 'r'
         <> metavar "[106] | INT"
         <> value 106
         <> help "Only search for stops with removeStops up to this\
                 \ amino acid position" )
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
    -- Get rid of carriages
    let contentsNoCarriages  = filter (/= '\r') $ unfilteredContents
    -- No newlines in sequence
    let contents              = joinSeq (removeN opts) contentsNoCarriages

    let genUnit               = aminoAcids opts
    let stopRange             = inputStopRange opts
    let codonMut              = inputCodonMut opts

    let cloneMap              = generateCloneMap contents
    let cloneMapNoStops       = if (removeStops opts)
                                    then removeStopsCloneMap
                                         genUnit stopRange cloneMap
                                    else cloneMap
    let cloneMapNoEmptyClones = removeEmptyClone cloneMapNoStops

    let outputString = if (removeGermlines opts)
                           then  printFastaNoGermline cloneMapNoEmptyClones
                           else  printFasta cloneMapNoEmptyClones

    writeFile (output opts) outputString

main :: IO ()
main = execParser opts >>= filterCLIPFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Filter CLIP fasta files in several optional ways"
     <> header "Filter CLIP Fasta, Gregory W. Schwartz" )
