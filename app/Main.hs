-- modify-fasta
-- By Gregory W. Schwartz

-- Takes a fasta file filters the fasta file in several ways.

{-# LANGUAGE BangPatterns #-}

-- Built-in
import Data.Maybe
import qualified Data.Map as M
import qualified System.IO as IO
import Control.Monad
import Control.Monad.Trans.Maybe

-- Cabal
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Options.Applicative
import Data.Fasta.Text
import Pipes
import qualified Pipes.Prelude as P
import qualified Pipes.Text as PT
import qualified Pipes.Text.IO as PT
import qualified Data.List.Split as Split
import Data.Semigroup ((<>))

-- Local
import Types
import Utility
import FilterCloneMap
import FilterFastaList
import FilterCloneList
import TransformFastaList
import TransformCloneList
import Print

-- Command line arguments
data Options = Options { input                    :: String
                       , aminoAcidsFlag           :: GeneticUnit
                       , legacyFlag               :: Bool
                       , clipFastaFlag            :: Bool
                       , convertToAminoAcidsFlag  :: Bool
                       , complementFlag           :: Bool
                       , reverseComplementFlag    :: Bool
                       , inputCodonTable          :: CodonTable
                       , inputFillIn              :: FillInValue
                       , inputStart               :: Maybe Int
                       , inputStop                :: Maybe Int
                       , inputMinSequenceMutation :: Maybe Int
                       , inputMutationCount       :: Maybe Int
                       , inputMutationPercent     :: Maybe Double
                       , addLengthFlag            :: Bool
                       , removeTheNsFlag          :: Bool
                       , removeGermlinesPreFlag   :: Bool
                       , removeHighlyMutatedFlag  :: Bool
                       , removeStopsFlag          :: Bool
                       , removeDuplicatesFlag     :: Bool
                       , removeOutOfFrameFlag     :: Bool
                       , removeUnknownNuc         :: Bool
                       , inputInFrame             :: Maybe Field
                       , inputOutFrame            :: Maybe Field
                       , trimFrame                :: Bool
                       , trackMutations           :: Maybe Field
                       , translateTrackMutations  :: Bool
                       , inputStopRange           :: Int
                       , inputCodonMut            :: CodonMut
                       , inputCodonMutType        :: String
                       , inputMutType             :: String
                       , inputChangeField         :: String
                       , inputCustomFilter        :: String
                       , customRemoveFlag         :: Bool
                       , customGermlineFlag       :: Bool
                       , inputInsertionPosition   :: Int
                       , inputReference           :: Maybe String
                       , inputGeneAlleleField     :: Int
                       , countFlag                :: Bool
                       , output                   :: String
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
      <*> option auto
          ( long "unit"
         <> short 'u'
         <> metavar "AminoAcid | Nucleotide"
         <> help "Whether these sequences are composed of\
                 \ amino acids (AminoAcid) or nucleotides (Nucleotide)" )
      <*> switch
          ( long "legacy"
         <> short 'L'
         <> help "Whether to use the legacy version with no pipes. Note: The\
                 \ legacy version will be unsupported and is greedy\
                 \ in terms of speed and memory. Use only if really needed.\
                 \ Features that are legacy only are noted in this\
                 \ documentation" )
      <*> switch
          ( long "clip-fasta"
         <> short 'A'
         <> help "Whether the input is a clip fasta file (has germline >>\
                 \ sequences)." )
      <*> switch
          ( long "convert-to-amino-acids"
         <> short 'C'
         <> help "Whether to convert the filtered sequences to amino acids\
                 \ in the output. Applied last, even after add length." )
      <*> switch
          ( long "complement"
         <> short 'b'
         <> help "Whether to find the complement of the sequence." )
      <*> switch
          ( long "reverse-complement"
         <> short 'B'
         <> help "Whether to to find the reverse complement of the sequence." )
      <*> option auto
          ( long "codon-table"
         <> metavar "[(CODON, AA)]"
         <> value []
         <> help "The codon table used for converting to amino acids. This\
                 \ is a list of codon AA pairs, so for example:\
                 \ [(\"AAT\", '@'), (\"ATC\", '#')] for AATATCGGG would result\
                 \ in @#G as defaults are used if the codon is not in\
                 \ the custom table." )
      <*> option auto
          ( long "fill-in"
         <> short 'F'
         <> metavar "(FIELD, START, 'CHARACTER')"
         <> value (-1, -1, 'X')
         <> help "Use the FIELD index (1 indexed, split by '|') in the\
                 \ header to fill in the unknown character CHARACTER with the\
                 \ respective character in that field. Use the START value\
                 \ to let the program know where the string in the field starts\
                 \ from in the sequence. For instance, a header of >H2O|HELLO\
                 \ for the sequence HIMANHELXODUDE and a value of fill-in of\
                 \ (2, 6, 'X') would change the sequence to HIMANHELLODUDE. If\
                 \ the string in the field has the unknown character as well,\
                 \ i.e. >H2O|HELXO, then the sequence is considered bad and is\
                 \ removed" )
      <*> optional ( option auto
          ( long "start"
         <> short 't'
         <> metavar "[ ] | INT"
         <> help "Remove everything before this position (1 indexed).\
                 \ Done first just after filtering." ) )
      <*> optional ( option auto
          ( long "stop"
         <> short 'p'
         <> metavar "[ ] | INT"
         <> help "Remove everything after this position (1 indexed).\
                 \ Done first just after filtering." ) )
      <*> optional ( option auto
          ( long "frequent-mutation-min"
         <> metavar "[ ] | INT"
         <> help "Minimum number of sequences required for a clone to be valid\
                 \ in the calculation of frequent mutations, replaces\
                 \ with gaps otherwise" ) )
      <*> optional ( option auto
          ( long "frequent-mutation-count"
         <> metavar "[ ] | INT"
         <> help "Only include codons containing a mutation present in this\
                 \ many sequences in the clone or more. 0 is all sequences.\
                 \ Converts the unincluded codons to gaps." ) )
      <*> optional ( option auto
          ( long "frequent-mutation-percent"
         <> metavar "[ ] | PERCENT"
         <> help "Only include codons containing a mutation present in this\
                 \ percentage of sequences in the clone or more.\
                 \ 0 is all sequences. Converts the unincluded codons\
                 \ to gaps." ) )
      <*> switch
          ( long "add-length"
         <> short 'l'
         <> help "Whether to append the length of the sequence to the end of\
                 \ the header, calculated after --convert-to-amino-acids\
                 \ if enabled" )
      <*> switch
          ( long "remove-N"
         <> short 'N'
         <> help "Whether to replace N or n in the sequence with a gap, '-'" )
      <*> switch
          ( long "remove-germlines"
         <> short 'g'
         <> help "Whether to remove germlines." )
      <*> switch
          ( long "remove-highly-mutated"
         <> short 'h'
         <> help "Whether to remove highly mutated clone sequences (a third\
                 \ of their sequence are different amino acids). Optionally\
                 \ uses codon-table." )
      <*> switch
          ( long "remove-stops"
         <> short 's'
         <> help "Whether to remove sequences with stop codons (\"*\").\
                 \ Optionally uses codon-table." )
      <*> switch
          ( long "legacy-remove-duplicates"
         <> short 'd'
         <> help "Whether to remove duplicate sequences. LEGACY ONLY" )
      <*> switch
          ( long "remove-out-of-frame"
         <> short 'O'
         <> help "Whether to remove sequences that are out of frame--if the\
                 \ sequences or number of gaps is not divisible by 3" )
      <*> switch
          ( long "remove-unknown-nucleotides"
         <> short 'Y'
         <> help "Convert unknown nucleotides (not ACGTN-.) to gaps (-)" )
      <*> optional ( option auto
          ( long "input-inframe"
         <> metavar "[ ] | FIELD"
         <> help "Represents the 1 indexed field split by '|'\
                 \ containing the inframe\
                 \ value (frames are 0, 1, or 2 like\
                 \ USCS definitions). For use with trim-frame." )
        )
      <*> optional ( option auto
          ( long "input-outframe"
         <> metavar "[ ] | FIELD"
         <> help "Represents the 1 indexed field split by '|'\
                 \ containing the outframe\
                 \ value (frames are 0, 1, or 2 like\
                 \ USCS definitions). For use with trim-frame." )
        )
      <*> switch
          ( long "trim-frame"
         <> short 'y'
         <> help "Trim each sequence to be in frame by remove extra nucleotides\
                 \ at the end. If input-inframe or input-outframe is\
                 \ specified, follow those rules instead." )
      <*> optional ( option auto
          ( long "track-mutations"
         <> metavar "[ ] | FIELD"
         <> help "Represents the 1 indexed field split by '|'\
                 \ containing the sequence to compare to. Takes the differences\
                 \ from the two sequences and puts them in the header." )
        )
      <*> switch
          ( long "translate-track-mutations"
         <> help "Whether to translate the field in track-mutations before\
                 \ finding differences." )
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
                 \ all codons). Converts the unincluded codon to gaps." )
      <*> strOption
          ( long "input-codon-mut-type"
         <> short 'T'
         <> metavar "[=]|>|<"
         <> value "="
         <> help "Only include codons with this many mutations (=)\
                 \ (or lesser (<) or greater (>), depending on\
                 \ input-codon-mut). Converts the unincluded codon to gaps.\
                 \ For use with input-codon-mut." )
      <*> strOption
          ( long "input-mut-type"
         <> short 'M'
         <> metavar "[All]|Silent|Replacement"
         <> value "All"
         <> help "Only include codons with all mutations (All),\
                 \ (or silent (Silent) or replacement (Replacement)).\
                 \ For use with input-codon-mut." )
      <*> strOption
          ( long "input-change-field"
         <> short 'e'
         <> metavar "(FIELD (Int), VALUE (String))"
         <> value ""
         <> help "Change a field to a match, so a regex \"ch.*_\" to field 2\
                 \ of \">abc|brie_cheese_dude\" would result in\
                 \ \">abc|cheese_\". Useful for getting specific properties\
                 \ from a field. Can take a list of format\
                 \ \"(Int, String)&&(Int, String)&& ...\" and so on. The String\
                 \ is in regex format (POSIX extended).\
                 \ The first in the tuple is the location of the field\
                 \ (1 indexed, split by '|')." )
      <*> strOption
          ( long "input-custom-filter"
         <> short 'f'
         <> metavar "(FIELD (Int), VALUE (String))"
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
          ( long "custom-remove"
         <> short 'm'
         <> help "Whether to remove the sequences containing the custom filter\
                 \ as opposed to remove the sequences that don't contain the\
                 \ filter" )
      <*> switch
          ( long "legacy-custom-germline"
         <> short 'G'
         <> help "Whether to apply the custom filter to germlines (>>)\
                 \ instead of sequences (>). LEGACY ONLY" )
      <*> option auto
          ( long "input-insertion-position"
         <> metavar "[1]|INT"
         <> value 1
         <> help "The field (1 indexed) of the where to insert the sequence\
                 \ into the reference from --input-insertion-reference.\
                 \ For instance, with a reference of \"ATT\", a position of 2,\
                 \ and a sequence of \"GC\", the result will be \"AGCTT\".\
                 \ Keeps the header of the insertion sequence. Requires\
                 \ --input-insertion-reference. No legacy." )
      <*> optional ( strOption
          ( long "input-reference"
         <> metavar "FILE"
         <> help "The file containing a reference sequence.\
                 \ Only the first sequence is used." )
          )
      <*> option auto
          ( long "legacy-gene-allele-field"
         <> short 'V'
         <> metavar "[1]|INT"
         <> value 1
         <> help "The field (1 indexed) of the gene allele name. LEGACY ONLY" )
      <*> switch
          ( long "legacy-count"
         <> short 'v'
         <> help "Do not save output, just count genes and alleles from\
                 \ the results. Requires gene-allele-field. LEGACY ONLY" )
      <*> strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> value ""
         <> help "The output fasta file" )

-- | Parse the argument fieldInt
fieldIntParser :: String -> [(Maybe Int, T.Text)]
fieldIntParser "" = []
fieldIntParser s  = map (\x -> (first x, second x)) . Split.splitOn "&&" $ s
  where
    first x
        | (head . Split.splitOn "," $ x) == "(" = Nothing
        | otherwise =
            Just (read (tail . head . Split.splitOn "," $ x) :: Int)
    second = T.pack . init . dropWhile (== ' ') . last . Split.splitOn ","

-- | Check for amino acid
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

    -- | Get a possible reference sequence.
    let getRefSeq :: MaybeT IO FastaSequence
        getRefSeq = do
            file   <- MaybeT . return $ inputReference opts
            hInRef <- lift $ IO.openFile file IO.ReadMode
            MaybeT . runEffect . P.head $ pipesFasta (PT.fromHandle hInRef)

    refSeq <- runMaybeT getRefSeq

    let genUnit        = aminoAcidsFlag opts
        stopRange      = inputStopRange opts
        customFilters  = fieldIntParser . inputCustomFilter $ opts
        changeFields   = fieldIntParser . inputChangeField $ opts
        codonMut       = inputCodonMut opts
        codonMutType   = T.pack . inputCodonMutType $ opts
        mutType        = T.pack . inputMutType $ opts
        insertPos      = inputInsertionPosition opts

        -- Remove out of frame sequences
        seqInFrame x = not ( removeOutOfFrameFlag opts
                          && (not . isAminoAcid $ genUnit)
                           )
                    || isInFrame x

        -- Start filtering out sequences
        -- Include only custom filter sequences
        customFilter x = null customFilters
                      || hasAllCustomFilters
                         (customRemoveFlag opts)
                         customFilters
                         x

        -- Remove clones with stops in the range
        noStops x = not (removeStopsFlag opts)
                 || hasNoStops genUnit (inputCodonTable opts) stopRange x

        -- Remove Ns from fasta list
        noNs = if removeTheNsFlag opts && (not . isAminoAcid $ genUnit)
                then removeN
                else id

        -- Change fasta headers with match
        changeHeader x = if not . null $ changeFields
                             then changeAllFields x changeFields
                             else x

        -- Get a specific region of the sequence
        cutSequence = case (inputStart opts, inputStop opts) of
                        (Nothing, Nothing) -> id
                        (start, stop)      -> getRegionSequence start stop

        -- Remove non standard nucleotides
        removeUnknown = if removeUnknownNuc opts
                            then removeUnknownNucs
                            else id

        -- Trim sequence
        trim fs =
            if trimFrame opts
                then trimFasta
                    genUnit
                    ((\x -> read . T.unpack . getField x '|' $ fs) <$> inputInFrame opts)
                    ((\x -> read . T.unpack . getField x '|' $ fs) <$> inputOutFrame opts)
                    fs
                else fs

        -- Put mutations in the header
        includeMutations fs =
            case trackMutations opts of
                Nothing  -> fs
                (Just x) ->
                    addMutationsHeader (translateTrackMutations opts) x fs

        -- Fill in bad characters at the requested section with possible
        -- replacements
        fillIn = case inputFillIn opts of
                    (-1, -1, 'X') -> id
                    (f, s, c)     -> fillInSequence f s c

        -- Insert a sequence into a reference sequence.
        insert fs = case (refSeq, insertPos) of
                        (Just r, p)  ->
                            insertSequence
                                (read . T.unpack . getField p '|' $ fs)
                                r
                                fs
                        _            -> fs

        -- Find the complement
        complement = if complementFlag opts
                        then compl
                        else id

        -- Convert to amino acids
        reverseComplement = if reverseComplementFlag opts
                                then revCompl
                                else id

        -- Convert to amino acids
        ntToaa = if convertToAminoAcidsFlag opts
                    then convertToAminoAcidsFastaSequence (inputCodonTable opts)
                    else id

        -- Include sequence length in header at the end
        includeLength = if addLengthFlag opts
                            then addLengthHeader
                            else id

        -- CLIP fasta specific filters and transformations

        -- Remove highly mutated sequences
        removeHighMutations = if removeHighlyMutatedFlag opts
                                then filterHighlyMutatedEntry
                                     genUnit
                                     (inputCodonTable opts)
                                else id

        -- Extract mutations to a certain degree
        getMutations = if codonMut > -1
                        then onlyMutations codonMut codonMutType mutType
                        else id

        -- Extract mutations to a certain degree
        getFrequentMutations = if isJust (inputMutationCount opts)
                               || isJust (inputMutationPercent opts)
                                   then frequentMutations
                                        (inputMinSequenceMutation opts)
                                        (inputMutationCount opts)
                                        (inputMutationPercent opts)
                                   else id

        -- Final order
        filterOrder x      = seqInFrame x && customFilter x && noStops x
        transformOrder     = includeLength
                           . includeMutations
                           . ntToaa
                           . reverseComplement
                           . complement
                           . changeHeader
                           . removeUnknown
                           . trim
                           . noNs
                           . fillIn
                           . insert
                           . cutSequence
        -- Specifically for germlines, as we don't want to change header or
        -- fill in the germline because that would make no sense in this
        -- case
        transformGermline  = includeLength
                           . includeMutations
                           . ntToaa
                           . reverseComplement
                           . complement
                           . removeUnknown
                           . trim
                           . noNs
                           . cutSequence
        -- Specifically for CLIP fasta files
        transformOrderCLIP = getFrequentMutations
                           . getMutations
                           . removeHighMutations

        executePrintFasta x = mappend (showFasta x) (T.pack "\n")
        executePrintCLIPFasta x =
            mappend
            (printCloneEntry (removeGermlinesPreFlag opts) x)
            (T.pack "\n")

    -- Execute pipes
    if clipFastaFlag opts
        then
            runEffect $ ( ( pipesCLIPFasta (PT.fromHandle hIn)
                        >-> P.map ( \(!germline, !fseqs) ->
                                    (germline, filter filterOrder fseqs)
                                  ) -- Filter sequences
                        >-> P.map transformOrderCLIP -- Transform specifically for CLIP fasta
                        >-> P.map ( \(!germline, !fseqs) ->
                                    ( transformGermline germline
                                    , map transformOrder fseqs
                                    )
                                  ) -- Transform sequences, only do some for germline
                        >-> P.map ( \(!germline, !fs)
                                 -> ( germline
                                    , filter (not . T.null . fastaSeq) fs
                                    )
                                  ) -- Remove empty sequences
                        >-> P.filter (not . null . snd) -- Remove empty clones
                        >-> P.map executePrintCLIPFasta ) -- Print the results
                         >> yield (T.pack "\n") )  -- want that newline at the end
                    >-> PT.toHandle hOut
        else
            runEffect $ ( ( pipesFasta (PT.fromHandle hIn)
                        >-> P.filter filterOrder -- Filter
                        >-> P.map transformOrder -- Transform
                        >-> P.filter (not . T.null . fastaSeq) -- Remove empty sequences
                        >-> P.map executePrintFasta ) -- Print the results
                         >> yield (T.pack "\n") )  -- want that newline at the end
                    >-> PT.toHandle hOut

    -- Finish up by closing file if written
    unless (null . input $ opts) (IO.hClose hIn)
    unless (null . output $ opts) (IO.hClose hOut)

-- Legacy function
modifyFastaCloneMap :: Options -> IO ()
modifyFastaCloneMap opts = do
    contents <- if null . input $ opts
                    then T.getContents
                    else T.readFile . input $ opts
    -- No redundant newlines in sequence
    let genUnit               = aminoAcidsFlag opts
        stopRange             = inputStopRange opts
        codonMut              = inputCodonMut opts
        codonMutType          = T.pack . inputCodonMutType $ opts
        mutType               = T.pack . inputMutType $ opts
        customFilters         = fieldIntParser $ inputCustomFilter opts
        removeGermlinesFlag   = if not . clipFastaFlag $ opts
                                    then True
                                    else removeGermlinesPreFlag opts

    -- Initiate CloneMap
        cloneMapFrames = if not . clipFastaFlag $ opts
                            then addFillerGermlines
                               . parsecFasta
                               $ contents
                            else parsecCLIPFasta contents

    -- Remove out of frame sequences
        cloneMapInFrame       = if (removeOutOfFrameFlag opts && ( not
                                                                 . isAminoAcid
                                                                 $ genUnit ) )
                                    then removeOutOfFrameSeqs cloneMapFrames
                                    else cloneMapFrames

    -- Remove Ns from CloneMap
        cloneMap              = if (removeTheNsFlag opts && ( not
                                                            . isAminoAcid
                                                            $ genUnit ) )
                                    then removeCLIPNs cloneMapInFrame
                                    else cloneMapInFrame

    -- Start filtering out sequences
    -- Include only custom filter sequences
        cloneMapCustom        = if not . null $ customFilters
                                    then removeAllCustomFilters
                                         (customGermlineFlag opts)
                                         (customRemoveFlag opts)
                                         cloneMap
                                         customFilters
                                    else cloneMap
    -- Remove clones with stops in the range
        (cloneMapNoStops, errorString) = if removeStopsFlag opts
                                            then removeStopsCloneMap
                                                 genUnit
                                                 (inputCodonTable opts)
                                                 stopRange
                                                 cloneMapCustom
                                            else (cloneMapCustom, Nothing)
    -- Output Error if necessary
    case errorString of
        Nothing -> return ()
        Just x  -> error x

    -- Remove duplicate sequences
    let cloneMapNoDuplicates  = if removeDuplicatesFlag opts
                                    then removeDuplicatesCloneMap
                                         cloneMapNoStops
                                    else cloneMapNoStops

    -- Remove clones that are highly mutated
        (cloneMapLowMutation, errorString2) = if (removeHighlyMutatedFlag opts)
                                              && ( not
                                                 . isAminoAcid
                                                 $ genUnit )
                                                 then filterHighlyMutated
                                                      genUnit
                                                      (inputCodonTable opts)
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
        (cloneMapAA, errorString3) = if (convertToAminoAcidsFlag opts)
                                     && (not . isAminoAcid $ genUnit)
                                      then convertToAminoAcidsCloneMap
                                           (inputCodonTable opts)
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
    case countFlag opts of
        True -> do
            -- Print results
            let outputText = printSequenceCount
                                (clipFastaFlag opts)
                                (inputGeneAlleleField opts)
                                cloneMapAA
            -- Print results to stdout
            T.putStrLn outputText
        False -> do
            -- Print results
            let outputText = if removeGermlinesFlag
                                then  printFastaNoGermline cloneMapAA
                                else  printFasta cloneMapAA

            -- Save results
            if null . output $ opts
                then T.putStrLn outputText
                else T.writeFile (output opts) outputText

modifyFasta :: Options -> IO ()
modifyFasta opts = if legacyFlag opts
                    then modifyFastaCloneMap opts
                    else modifyFastaList opts

main :: IO ()
main = execParser opts >>= modifyFasta
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Modify fasta (and CLIP) files in several optional ways.\
                 \ Order of transformation goes: seqInFrame -> customFilter\
                 \ -> noStops -> removeHighMutations -> getMutations ->\
                 \ getFrequentMutations -> cutSequence -> insert -> fillIn\
                 \ -> noNs -> changeHeader -> complement -> reverseComplement\
                 \ -> ntToaa -> includeLength,\
                 \ so if you require a different\
                 \ order (which can change results dramatically), then do\
                 \ so one at a time through the wonderful world of piping."
     <> header "modify-fasta, Gregory W. Schwartz" )
