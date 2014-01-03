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
data Options = Options { input               :: String
                       , aminoAcids          :: GeneticUnit
                       , normalFasta         :: Bool
                       , convertToAminoAcids :: Bool
                       , removeN             :: Bool
                       , removeGermlines     :: Bool
                       , removeStops         :: Bool
                       , inputStopRange      :: Int
                       , inputCodonMut       :: CodonMut
                       , inputCustomFilter   :: String
                       , inputCustomField    :: Int
                       , infixCustomFilter   :: Bool
                       , customGermline      :: Bool
                       , customRemove        :: Bool
                       , output              :: String
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
          ( long "normalFasta"
         <> short 'A'
         <> help "Whether the input is a normal fasta file (no germline >>\
                 \ sequences" )
      <*> switch
          ( long "convertToAminoAcids"
         <> short 'C'
         <> help "Whether to convert the filtered sequences to amino acids\
                 \ in the output" )
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
         <> metavar "[106]|INT"
         <> value 106
         <> help "Only search for stops with removeStops up to this\
                 \ amino acid position" )
      <*> option
          ( long "inputCodonMut"
         <> short 'c'
         <> metavar "[0]|1|2|3"
         <> value 0
         <> help "Only include codons with this many mutations or less\
                 \ (0 is the same as include all codons). Converts\
                 \ the codon to gaps" )
      <*> strOption
          ( long "inputCustomFilter"
         <> short 'f'
         <> metavar "FIELD_VALUE (String)"
         <> value ""
         <> help "A custom field value to keep, discarding all others" )
      <*> option
          ( long "inputCustomField"
         <> short 'F'
         <> metavar "FIELD_LOCATION (Int)"
         <> value 0
         <> help "The location of the field in inputCustomFilter\
                 \ (split by |, 0 means search the whole header)" )
      <*> switch
          ( long "infixCustomFilter"
         <> short 'I'
         <> help "Whether to find the custom filter in the field\
                 \ as an infix (some part of the custom field\
                 \ matches the header) as opposed to an\
                 \ exact match (the entire field must be the custom field)" )
      <*> switch
          ( long "customGermline"
         <> short 'G'
         <> help "Whether to apply the custom filter to germlines (>>)\
                 \ instead of sequences (>)" )
      <*> switch
          ( long "customRemove"
         <> short 'm'
         <> help "Whether to remove the sequences containing the custom filter\
                 \ as opposed to remove the sequences that don't contain the\
                 \ filter" )
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
    let contentsNoCarriages   = filter (/= '\r') $ unfilteredContents
    -- If normalFasta, insert filler germlines
    let contentsCLIP          = if (normalFasta opts)
                                    then addFillerGermlines contentsNoCarriages
                                    else contentsNoCarriages
    -- No redundant newlines in sequence
    let contents              = joinSeq (removeN opts) contentsCLIP

    let genUnit               = aminoAcids opts
    let stopRange             = inputStopRange opts
    let codonMut              = inputCodonMut opts
    let customFilter          = inputCustomFilter opts
    let customField           = inputCustomField opts
    let removeGermlinesFlag   = if (normalFasta opts)
                                    then True
                                    else (removeGermlines opts)

    -- Initiate CloneMap
    let cloneMap              = generateCloneMap contents

    -- Start filtering out sequences
    -- Include only custom filter sequences
    let cloneMapCustom        = if (not . null $ customFilter)
                                    then removeCustomFilter
                                         (customGermline opts)
                                         (customRemove opts)
                                         (infixCustomFilter opts)
                                         customField
                                         customFilter
                                         cloneMap
                                    else cloneMap
    -- Remove clones with stops in the range
    let cloneMapNoStops       = if (removeStops opts)
                                    then removeStopsCloneMap
                                         genUnit stopRange cloneMapCustom
                                    else cloneMapCustom
    -- Remove codons with codons with a certain number of mutations
    let cloneMapNoCodonMut    = if (codonMut > 0)
                                    then removeCodonMutCount
                                         codonMut cloneMapNoStops
                                    else cloneMapNoStops
    -- Remove empty clones
    let cloneMapNoEmptyClones = removeEmptyClone cloneMapNoCodonMut

    -- Convert sequences to amino acids
    let cloneMapAA            = if (convertToAminoAcids opts)
                                    then convertToAminoAcidsCloneMap
                                         cloneMapNoEmptyClones
                                    else cloneMapNoEmptyClones

    -- Print results
    let outputString = if removeGermlinesFlag
                           then  printFastaNoGermline cloneMapAA
                           else  printFasta cloneMapAA

    -- Save results
    writeFile (output opts) outputString

main :: IO ()
main = execParser opts >>= filterCLIPFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Filter CLIP fasta files in several optional ways"
     <> header "Filter CLIP Fasta, Gregory W. Schwartz" )
