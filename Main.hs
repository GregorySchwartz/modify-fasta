-- Modify Fasta
-- By G.W. Schwartz

-- Takes a CLIP fasta file or fasta file and filters the fasta file in
-- several optional ways.

-- Cabal
import Options.Applicative

-- Local
import Types
import FastaParse
import FilterCloneMap
import Print
import Control.Monad.State
import qualified Data.List.Split as Split

-- Command line arguments
data Options = Options { input               :: String
                       , aminoAcids          :: GeneticUnit
                       , normalFasta         :: Bool
                       , convertToAminoAcids :: Bool
                       , removeN             :: Bool
                       , removeGermlines     :: Bool
                       , removeHighlyMutated :: Bool
                       , removeStops         :: Bool
                       , inputStopRange      :: Int
                       , inputCodonMut       :: CodonMut
                       , inputCodonMutType   :: String
                       , inputMutType        :: String
                       , inputCustomFilter   :: String
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
                 \ sequences)" )
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
          ( long "removeHighlyMutated"
         <> short 'h'
         <> help "Whether to remove highly mutated clone sequences (a third\
                 \ of their sequence are different amino acids)" )
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
         <> metavar "[-1]|0|1|2|3"
         <> value (-1)
         <> help "Only include codons with this many mutations or less or more,\
                 \ depending on inputCodonMutType (-1 is the same as include\
                 \ all codons). Converts the codon to gaps" )
      <*> strOption
          ( long "inputCodonMutType"
         <> short 'T'
         <> metavar "[=]|>|<"
         <> value "="
         <> help "Only include codons with this many mutations (=)\
                 \ (or lesser (<) or greater (>), depending on\
                 \ inputCodonMut). Converts the codon to gaps" )
      <*> strOption
          ( long "inputMutType"
         <> short 'M'
         <> metavar "[All]|Silent|Replacement"
         <> value "All"
         <> help "Only include codons with this all mutations (All),\
                 \ (or silent (Silent) or replacement (Replacement))" )
      <*> strOption
          ( long "inputCustomFilter"
         <> short 'f'
         <> metavar "((FIELD_LOCATION (Int), FIELD_VALUE (String))"
         <> value ""
         <> help "A custom filter. Can take a list of format\
                 \ \"(Int, String)&&(Int, String)&& ...\" and so on.\
                 \ The first in the tuple is the location of the field\
                 \ (1 indexed, split by '|'). If you want to apply to\
                 \ the entire header, either have the location as 0 or\
                 \ exclude the location altogether (, Day 3|IGHV3) for instance\
                 \ will match if the entire header is '>Day 3|IGHV3'.\
                 \ This list will be filtered one at a time, so you cannot\
                 \ get multiple filters, but you can remove multiple filters." )
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

customFiltersIntParser :: String -> [(Maybe Int, String)]
customFiltersIntParser "" = []
customFiltersIntParser s = map (\x -> (first x, second x))
                         . Split.splitOn "&&"
                         $ s
  where
    first x
        | (head . Split.splitOn "," $ x) == "(" = Nothing
        | otherwise =
            Just (read (tail . head . Split.splitOn "," $ x) :: Int)
    second = init . dropWhile (== ' ') . last . Split.splitOn ","

modifyFasta :: Options -> IO ()
modifyFasta opts = do
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
    let codonMutType          = inputCodonMutType opts
    let mutType               = inputMutType opts
    let customFilters         = customFiltersIntParser $ inputCustomFilter opts
    let removeGermlinesFlag   = if (normalFasta opts)
                                    then True
                                    else (removeGermlines opts)

    -- Initiate CloneMap
    let cloneMap              = generateCloneMap contents

    -- Start filtering out sequences
    -- Include only custom filter sequences
    let cloneMapCustom        = if (not . null $ customFilters)
                                    then snd
                                       . runState ( removeAllCustomFilters
                                                    (customGermline opts)
                                                    (customRemove opts)
                                                    (infixCustomFilter opts)
                                                    customFilters )
                                       $ cloneMap
                                    else cloneMap
    -- Remove clones with stops in the range
    let cloneMapNoStops       = if (removeStops opts)
                                    then removeStopsCloneMap
                                         genUnit stopRange cloneMapCustom
                                    else cloneMapCustom
    -- Remove clones that are highly mutated
    let cloneMapLowMutation   = if (removeHighlyMutated opts)
                                    then filterHighlyMutated
                                         genUnit cloneMapCustom
                                    else cloneMapCustom
    -- Remove codons with codons with a certain number of mutations
    let cloneMapNoCodonMut    = if (codonMut > -1)
                                    then removeCodonMutCount codonMut
                                                             codonMutType
                                                             mutType
                                                             cloneMapLowMutation
                                    else cloneMapLowMutation
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
main = execParser opts >>= modifyFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Modify fasta (and CLIP) files in several optional ways"
     <> header "Modify Fasta, Gregory W. Schwartz" )
