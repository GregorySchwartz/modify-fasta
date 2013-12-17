-- Remove Germline
-- By G.W. Schwartz

-- Takes a fasta CLIP file and converts the file into a normal fasta file
-- by removing the germline sequences (and thus information about clones is
-- lost).

-- Cabal
import Options.Applicative

-- Local
import Types
import FastaParse
import Print

-- Command line arguments
data Options = Options { input  :: String
                       , output :: String
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
      <*> option
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> value "output.fasta"
         <> help "The output fasta file" )

removeGermline :: Options -> IO ()
removeGermline opts = do
    contents <- readFile . input $ opts

    let unfilteredFastaList = fastaParser contents
    let fastaList           = filterFasta unfilteredFastaList

    writeFile (output opts) . printFasta $ fastaList

main :: IO ()
main = execParser opts >>= removeGermline
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Remove the germline sequences from a CLIP fasta file"
     <> header "Remove Germline, Gregory W. Schwartz" )
