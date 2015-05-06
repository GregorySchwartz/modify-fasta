-- Modify Fasta
-- By G.W. Schwartz

-- Takes a CLIP fasta file or fasta file and filters the fasta file in
-- several optional ways.

-- Built-in
import qualified Data.Map as M
import qualified System.IO as IO
import Control.Monad

-- Cabal
import Options.Applicative
import Data.Fasta.String
import Pipes
import qualified Pipes.Prelude as P
import qualified Data.List.Split as Split

-- Local
import Types
import Utility
import FilterCloneMap
import FilterFastaList
import Print

-- Command line arguments
data Options = Options { input               :: String
                       , aminoAcids          :: String
                       , legacy              :: Bool
                       , clipFasta           :: Bool
                       , convertToAminoAcids :: Bool
                       , removeTheNs         :: Bool
                       , removeGermlines     :: Bool
                       , removeHighlyMutated :: Bool
                       , removeStops         :: Bool
                       , removeDuplicates    :: Bool
                       , removeOutOfFrame    :: Bool
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
      <*> strOption
          ( long "unit"
         <> short 'u'
         <> metavar "AminoAcid|Nucleotide"
         <> help "Whether these sequences are composed of\
                 \ amino acids (AminoAcid) or nucleotides (Nucleotide)" )
      <*> switch
          ( long "legacy"
         <> short 'L'
         <> help "Whether to use the legacy version with no pipes. Note: The\
                 \ legacy version supports more features but is greedy\
                 \ in terms of speed and memory. Use only if really needed.\
                 \ Features that are legacy only are noted in this\
                 \ documentation" )
      <*> switch
          ( long "clip-fasta"
         <> short 'A'
         <> help "Whether the input is a clip fasta file (has germline >>\
                 \ sequences). LEGACY ONLY" )
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
         <> help "Whether to remove germlines. LEGACY ONLY" )
      <*> switch
          ( long "remove-highly-mutated"
         <> short 'h'
         <> help "Whether to remove highly mutated clone sequences (a third\
                 \ of their sequence are different amino acids). LEGACY ONLY" )
      <*> switch
          ( long "remove-stops"
         <> short 's'
         <> help "Whether to remove sequences with stop codons" )
      <*> switch
          ( long "remove-duplicates"
         <> short 'd'
         <> help "Whether to remove duplicate sequences. LEGACY ONLY" )
      <*> switch
          ( long "remove-out-of-frame"
         <> short 'O'
         <> help "Whether to remove sequences that are out of frame--if the\
                 \ sequences or number of gaps is not divisible by 3" )
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
                 \ all codons). Converts the codon to gaps. LEGACY ONLY" )
      <*> strOption
          ( long "input-codon-mut-type"
         <> short 'T'
         <> metavar "[=]|>|<"
         <> value "="
         <> help "Only include codons with this many mutations (=)\
                 \ (or lesser (<) or greater (>), depending on\
                 \ input-codon-mut). Converts the codon to gaps. LEGACY ONLY" )
      <*> strOption
          ( long "input-mut-type"
         <> short 'M'
         <> metavar "[All]|Silent|Replacement"
         <> value "All"
         <> help "Only include codons with this all mutations (All),\
                 \ (or silent (Silent) or replacement (Replacement)). LEGACY\
                 \ ONLY" )
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
                 \ instead of sequences (>). LEGACY ONLY" )
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
         <> help "The field (1 indexed) of the gene allele name. LEGACY ONLY" )
      <*> switch
          ( long "count"
         <> short 'v'
         <> help "Do not save output, just count genes and alleles from\
                 \ the results. Requires gene-allele-field. LEGACY ONLY" )
      <*> strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> value ""
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

modifyFastaList :: Options -> IO ()
modifyFastaList opts = do
    hIn  <- if null . input $ opts
                then return IO.stdin
                else IO.openFile (input opts) IO.ReadMode
    hOut <- if null . output $ opts
                then return IO.stdout
                else IO.openFile (output opts) IO.WriteMode
    let genUnit               = read . aminoAcids $ opts
        stopRange             = inputStopRange opts
        customFilters         = customFiltersIntParser
                              . inputCustomFilter
                              $ opts

    -- Remove out of frame sequences
        seqInFrame x = if removeOutOfFrame opts && (not . isAminoAcid $ genUnit)
                        then isInFrame x
                        else True

    -- Remove Ns from CloneMap
        noNs x = if removeTheNs opts && (not . isAminoAcid $ genUnit)
                    then removeN x
                    else x
    -- Start filtering out sequences
    -- Include only custom filter sequences
        customFilter x = if not . null $ customFilters
                            then hasAllCustomFilters
                                 (customRemove opts)
                                 customFilters
                                 x
                            else True

    -- Remove clones with stops in the range
        noStops x = if removeStops opts
                        then hasNoStops
                             genUnit
                             stopRange
                             x
                        else True

    -- Filter
    runEffect $ P.fromHandle hIn
            >-> pipesFasta hIn
            >-> P.filter (\x -> seqInFrame x && customFilter x && noStops x)
            >-> P.map ( showFasta
                      . (\x -> if convertToAminoAcids opts then convertToAminoAcidsFastaSequence x else x)
                      . noNs )
            >-> P.toHandle hOut

    -- Finish up by closing file if written
    unless (null . output $ opts) (IO.hClose hOut)

modifyFastaCloneMap :: Options -> IO ()
modifyFastaCloneMap opts = do
    contents <- if null . input $ opts
                    then getContents
                    else readFile . input $ opts
    -- No redundant newlines in sequence
    let genUnit               = read . aminoAcids $ opts
        stopRange             = inputStopRange opts
        codonMut              = inputCodonMut opts
        codonMutType          = inputCodonMutType opts
        mutType               = inputMutType opts
        customFilters         = customFiltersIntParser $ inputCustomFilter opts
        removeGermlinesFlag   = if not . clipFasta $ opts
                                    then True
                                    else removeGermlines opts

    -- Initiate CloneMap
        cloneMapFrames = if not . clipFasta $ opts
                            then addFillerGermlines
                               . parseFasta
                               $ contents
                            else parseCLIPFasta contents

    -- Remove out of frame sequences
        cloneMapInFrame       = if (removeOutOfFrame opts && ( not
                                                             . isAminoAcid
                                                             $ genUnit ) )
                                    then removeOutOfFrameSeqs cloneMapFrames
                                    else cloneMapFrames

    -- Remove Ns from CloneMap
        cloneMap              = if (removeTheNs opts && ( not
                                                        . isAminoAcid
                                                        $ genUnit ) )
                                    then removeCLIPNs cloneMapInFrame
                                    else cloneMapInFrame

    -- Start filtering out sequences
    -- Include only custom filter sequences
        cloneMapCustom        = if not . null $ customFilters
                                    then removeAllCustomFilters
                                         (customGermline opts)
                                         (customRemove opts)
                                         cloneMap
                                         customFilters
                                    else cloneMap
    -- Remove clones with stops in the range
        (cloneMapNoStops, errorString) = if removeStops opts
                                            then removeStopsCloneMap
                                                 genUnit
                                                 stopRange
                                                 cloneMapCustom
                                            else (cloneMapCustom, Nothing)
    -- Output Error if necessary
    case errorString of
        Nothing -> return ()
        Just x  -> error x

    -- Remove duplicate sequences
    let cloneMapNoDuplicates  = if removeDuplicates opts
                                    then removeDuplicatesCloneMap
                                         cloneMapNoStops
                                    else cloneMapNoStops

    -- Remove clones that are highly mutated
        (cloneMapLowMutation, errorString2) = if (removeHighlyMutated opts)
                                              && ( not
                                                 . isAminoAcid
                                                 $ genUnit )
                                                 then filterHighlyMutated
                                                      genUnit
                                                      cloneMapNoDuplicates
                                                 else ( cloneMapNoDuplicates
                                                      , Nothing )

    -- Output Error if necessary
    case errorString2 of
        Nothing -> return ()
        Just x  -> error x

    -- Remove codons with codons with a certain number of mutations
    let cloneMapNoCodonMut    = if codonMut > -1
                                    then removeCodonMutCount codonMut
                                                             codonMutType
                                                             mutType
                                                             cloneMapLowMutation
                                    else cloneMapLowMutation
    -- Remove empty clones
        cloneMapNoEmptyClones = removeEmptyClone cloneMapNoCodonMut

    -- Convert sequences to amino acids
        (cloneMapAA, errorString3) = if (convertToAminoAcids opts)
                                     && (not . isAminoAcid $ genUnit)
                                      then convertToAminoAcidsCloneMap
                                           cloneMapNoEmptyClones
                                      else (cloneMapNoEmptyClones, Nothing)

    -- Output Error if necessary
    case errorString3 of
        Nothing -> return ()
        Just x  -> error x

    -- Break if there are no sequences to output
    case M.null cloneMapAA of
        True -> error "No sequences left! Nothing written."
        False -> return ()

    -- What to do with results
    case count opts of
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
            if null . output $ opts
                then putStrLn outputString
                else writeFile (output opts) outputString

modifyFasta :: Options -> IO ()
modifyFasta opts = if legacy opts
                    then modifyFastaCloneMap opts
                    else modifyFastaList opts

main :: IO ()
main = execParser opts >>= modifyFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Modify fasta (and CLIP) files in several optional ways"
     <> header "Modify Fasta, Gregory W. Schwartz" )
