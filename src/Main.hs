-- Modify Fasta
-- By G.W. Schwartz

-- Takes a CLIP fasta file or fasta file and filters the fasta file in
-- several optional ways.

-- Cabal
import Options.Applicative
import Data.Fasta.String.Parse

-- Local
import Types
import Utility
import FilterCloneMap
import Print
import qualified Data.List.Split as Split

-- Command line arguments
data Options = Options { input               :: String
                       , aminoAcids          :: GeneticUnit
                       , clipFasta           :: Bool
                       , convertToAminoAcids :: Bool
                       , removeN             :: Bool
                       , removeGermlines     :: Bool
                       , removeHighlyMutated :: Bool
                       , removeStops         :: Bool
                       , removeDuplicates    :: Bool
                       , inputStopRange      :: Int
                       , inputCodonMut       :: CodonMut
                       , inputCodonMutType   :: String
                       , inputMutType        :: String
                       , inputCustomFilter   :: String
                       , customGermline      :: Bool
                       , customRemove        :: Bool
                       , geneAlleleField     :: Int
                       , count               :: Bool
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
         <> help "The input fasta file or CLIP fasta file" )
      <*> flag Nucleotide AminoAcid
          ( long "amino-acids"
         <> short 'a'
         <> help "Whether these sequences are composed of\
                 \ amino acids (AminoAcid) or nucleotides (Nucleotide)" )
      <*> switch
          ( long "clip-fasta"
         <> short 'A'
         <> help "Whether the input is a clip fasta file (has germline >>\
                 \ sequences)" )
      <*> switch
          ( long "convert-to-amino-acids"
         <> short 'C'
         <> help "Whether to convert the filtered sequences to amino acids\
                 \ in the output" )
      <*> switch
          ( long "remove-N"
         <> short 'N'
         <> help "Whether to remove N or n in the sequence" )
      <*> switch
          ( long "remove-germlines"
         <> short 'g'
         <> help "Whether to remove germlines" )
      <*> switch
          ( long "remove-highly-mutated"
         <> short 'h'
         <> help "Whether to remove highly mutated clone sequences (a third\
                 \ of their sequence are different amino acids)" )
      <*> switch
          ( long "remove-stops"
         <> short 's'
         <> help "Whether to remove sequences with stop codons" )
      <*> switch
          ( long "remove-duplicates"
         <> short 'd'
         <> help "Whether to remove duplicate sequences" )
      <*> option auto
          ( long "input-stop-range"
         <> short 'r'
         <> metavar "[106]|INT"
         <> value 106
         <> help "Only search for stops with remove-stops up to this\
                 \ amino acid position" )
      <*> option auto
          ( long "input-codon-mut"
         <> short 'c'
         <> metavar "[-1]|0|1|2|3"
         <> value (-1)
         <> help "Only include codons with this many mutations or less or more,\
                 \ depending on input-codon-mut-type (-1 is the same as include\
                 \ all codons). Converts the codon to gaps" )
      <*> strOption
          ( long "input-codon-mut-type"
         <> short 'T'
         <> metavar "[=]|>|<"
         <> value "="
         <> help "Only include codons with this many mutations (=)\
                 \ (or lesser (<) or greater (>), depending on\
                 \ input-codon-mut). Converts the codon to gaps" )
      <*> strOption
          ( long "input-mut-type"
         <> short 'M'
         <> metavar "[All]|Silent|Replacement"
         <> value "All"
         <> help "Only include codons with this all mutations (All),\
                 \ (or silent (Silent) or replacement (Replacement))" )
      <*> strOption
          ( long "input-custom-filter"
         <> short 'f'
         <> metavar "((FIELD_LOCATION (Int), FIELD_VALUE (String))"
         <> value ""
         <> help "A custom filter. Can take a list of format\
                 \ \"(Int, String)&&(Int, String)&& ...\" and so on. The String\
                 \ is in regex format (POSIX extended), so if the entire\
                 \ string 'a' is in the whole field, then you need to input\
                 \ '^a$' for the beginning to the end!\
                 \ The first in the tuple is the location of the field\
                 \ (1 indexed, split by '|'). If you want to apply to\
                 \ the entire header, either have the location as 0 or\
                 \ exclude the location altogether (, Day 3|IGHV3) for instance\
                 \ will match if the entire header is '>Day 3|IGHV3'.\
                 \ This list will be filtered one at a time, so you cannot\
                 \ get multiple filters, but you can remove multiple filters." )
      <*> switch
          ( long "custom-germline"
         <> short 'G'
         <> help "Whether to apply the custom filter to germlines (>>)\
                 \ instead of sequences (>)" )
      <*> switch
          ( long "custom-remove"
         <> short 'm'
         <> help "Whether to remove the sequences containing the custom filter\
                 \ as opposed to remove the sequences that don't contain the\
                 \ filter" )
      <*> option auto
          ( long "gene-allele-field"
         <> short 'V'
         <> metavar "[1]|INT"
         <> value 1
         <> help "The field (1 indexed) of the gene allele name" )
      <*> switch
          ( long "count"
         <> short 'v'
         <> help "Do not save output, just count genes and alleles from\
                 \ the results. Requires gene-allele-field" )
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

isAminoAcid :: GeneticUnit -> Bool
isAminoAcid AminoAcid = True
isAminoAcid _         = False

modifyFasta :: Options -> IO ()
modifyFasta opts = do
    -- Checking if user knows what he or she is doing
    case (aminoAcids opts) of
        AminoAcid -> putStrLn "Amino acid sequence input, nucleotide options if\
                              \ applicable will be ignored"
        _         -> putStrLn "Nucleotide sequence input, amino acid options if\
                              \ applicable will be ignored"

    contents <- readFile . input $ opts
    -- No redundant newlines in sequence
    let genUnit               = aminoAcids opts
        stopRange             = inputStopRange opts
        codonMut              = inputCodonMut opts
        codonMutType          = inputCodonMutType opts
        mutType               = inputMutType opts
        customFilters         = customFiltersIntParser $ inputCustomFilter opts
        removeGermlinesFlag   = if (not . clipFasta $ opts)
                                    then True
                                    else (removeGermlines opts)

    -- Initiate CloneMap
        cloneMapNs            = if (not . clipFasta $ opts)
                                    then addFillerGermlines
                                       . parseFasta
                                       $ contents
                                    else parseCLIPFasta contents

    -- Remove Ns from CloneMap
        cloneMap              = if (removeN opts)
                                    then removeCLIPNs cloneMapNs
                                    else cloneMapNs

    -- Start filtering out sequences
    -- Include only custom filter sequences
        cloneMapCustom        = if (not . null $ customFilters)
                                    then removeAllCustomFilters
                                         (customGermline opts)
                                         (customRemove opts)
                                         cloneMap
                                         customFilters
                                    else cloneMap
    -- Remove clones with stops in the range
        (cloneMapNoStops, errorString) = if (removeStops opts)
                                            then removeStopsCloneMap
                                                 genUnit
                                                 stopRange
                                                 cloneMapCustom
                                            else (cloneMapCustom, Nothing)
    -- Output Error if necessary
    case errorString of
        Nothing -> return ()
        Just x  -> putStrLn x

    -- Remove duplicate sequences
    let cloneMapNoDuplicates  = if (removeDuplicates opts)
                                    then removeDuplicatesCloneMap
                                         cloneMapNoStops
                                    else cloneMapNoStops

    -- Remove clones that are highly mutated
        (cloneMapLowMutation, errorString2) = if (removeHighlyMutated opts)
                                              && ( not
                                                 . isAminoAcid
                                                 . aminoAcids
                                                 $ opts )
                                                 then filterHighlyMutated
                                                      genUnit
                                                      cloneMapNoDuplicates
                                                 else ( cloneMapNoDuplicates
                                                      , Nothing )

    -- Output Error if necessary
    case errorString2 of
        Nothing -> return ()
        Just x  -> putStrLn x

    -- Remove codons with codons with a certain number of mutations
    let cloneMapNoCodonMut    = if (codonMut > -1)
                                    then removeCodonMutCount codonMut
                                                             codonMutType
                                                             mutType
                                                             cloneMapLowMutation
                                    else cloneMapLowMutation
    -- Remove empty clones
        cloneMapNoEmptyClones = removeEmptyClone cloneMapNoCodonMut

    -- Convert sequences to amino acids
        (cloneMapAA, errorString3) = if (convertToAminoAcids opts)
                                     && (not . isAminoAcid . aminoAcids $ opts)
                                      then convertToAminoAcidsCloneMap
                                           cloneMapNoEmptyClones
                                      else (cloneMapNoEmptyClones, Nothing)

    -- Output Error if necessary
    case errorString3 of
        Nothing -> return ()
        Just x  -> putStrLn x

    -- What to do with results
    case (count opts) of
        True -> do
            -- Print results
            let outputString = printSequenceCount
                               (clipFasta opts)
                               (geneAlleleField opts)
                               cloneMapAA
            -- Print results to stdout
            putStrLn outputString
        False -> do
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
