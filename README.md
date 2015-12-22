modify-fasta
============

**Gregory W. Schwartz**

This program will take a fasta file, do certain filterings and transformations
as specified by the user, and return the new fasta file.

## Installation
```
stack install --resolver nightly
```

## Examples

The recommended usage for this program is to pipe in a fasta stream and make
only one transformation at a time (not necessary, but the ordering can get
confusing so to minimize assumptions about the ordering it's best to do it in
the order you request). For instance, to translate the sequences:

`cat input.fasta | modify-fasta -u Nucleotides --convert-to-amino-acids >
output.fasta`

To trim trailing nucleotides:

`cat input.fasta | modify-fasta -u Nucleotides --trim > output.fasta`

To only return `MOUSE` entries (where the entry header looks like
"acc15|inframe|3|MOUSE|False") and cut off the amino acids after position 20:

`cat input.fasta | modify-fasta -u AminoAcid --input-custom-filter "(4,
^MOUSE$)" | modify-fasta -u AminoAcid --stop 20 > output.fasta`

## Usage

```
modify-fasta, Gregory W. Schwartz

Usage: modify-fasta [-i|--input FILE] (-u|--unit AminoAcid|Nucleotide)
                    [-L|--legacy] [-A|--clip-fasta]
                    [-C|--convert-to-amino-acids]
                    [-F|--fill-in (FIELD, START, 'CHARACTER')]
                    [-t|--start [ ] | INT] [-p|--stop [ ] | INT]
                    [--frequent-mutation-min [ ] | INT]
                    [--frequent-mutation-count [ ] | INT]
                    [--frequent-mutation-percent [ ] | PERCENT]
                    [-l|--add-length] [-N|--remove-N] [-g|--remove-germlines]
                    [-h|--remove-highly-mutated] [-s|--remove-stops]
                    [-d|--legacy-remove-duplicates] [-O|--remove-out-of-frame]
                    [-Y|--remove-unknown-nucleotides]
                    [--input-inframe [ ] | FIELD] [--input-outframe [ ] | FIELD]
                    [-y|--trim-frame] [-r|--input-stop-range [106]|INT]
                    [-c|--input-codon-mut [-1]|0|1|2|3]
                    [-T|--input-codon-mut-type [=]|>|<]
                    [-M|--input-mut-type [All]|Silent|Replacement]
                    [-e|--input-change-field ((FIELD (Int), VALUE (String))]
                    [-f|--input-custom-filter ((FIELD (Int), VALUE (String))]
                    [-G|--legacy-custom-germline] [-m|--custom-remove]
                    [-V|--legacy-gene-allele-field [1]|INT] [-v|--legacy-count]
                    [-o|--output FILE]
  Modify fasta (and CLIP) files in several optional ways. Order of
  transformation goes: seqInFrame -> customFilter -> noStops ->
  removeHighMutations -> getMutations -> getFrequentMutations -> cutSequence ->
  fillIn -> noNs -> changeHeader -> ntToaa -> includeLength, so if you require a
  different order (which can change results dramatically), then do so one at a
  time through the wonderful world of piping.

Available options:
  -h,--help                Show this help text
  -i,--input FILE          The input fasta file or CLIP fasta file
  -u,--unit AminoAcid|Nucleotide
                           Whether these sequences are composed of amino acids
                           (AminoAcid) or nucleotides (Nucleotide)
  -L,--legacy              Whether to use the legacy version with no pipes.
                           Note: The legacy version supports more features but
                           is greedy in terms of speed and memory. Use only if
                           really needed. Features that are legacy only are
                           noted in this documentation
  -A,--clip-fasta          Whether the input is a clip fasta file (has germline
                           >> sequences).
  -C,--convert-to-amino-acids
                           Whether to convert the filtered sequences to amino
                           acids in the output. Applied last, even after add
                           length.
  -F,--fill-in (FIELD, START, 'CHARACTER')
                           Use the FIELD index (1 indexed, split by '|') in the
                           header to fill in the unknown character CHARACTER
                           with the respective character in that field. Use the
                           START value to let the program know where the string
                           in the field starts from in the sequence. For
                           instance, a header of >H2O|HELLO for the sequence
                           HIMANHELXODUDE and a value of fill-in of (2, 6, 'X')
                           would change the sequence to HIMANHELLODUDE. If the
                           string in the field has the unknown character as
                           well, i.e. >H2O|HELXO, then the sequence is
                           considered bad and is removed
  -t,--start [ ] | INT     Remove everything before this position (1 indexed).
                           Done first just after filtering.
  -p,--stop [ ] | INT      Remove everything after this position (1 indexed).
                           Done first just after filtering.
  --frequent-mutation-min [ ] | INT
                           Minimum number of sequences required for a clone to
                           be valid in the calculation of frequent mutations,
                           replaces with gaps otherwise
  --frequent-mutation-count [ ] | INT
                           Only include codons containing a mutation present in
                           this many sequences in the clone or more. 0 is all
                           sequences. Converts the unincluded codons to gaps.
  --frequent-mutation-percent [ ] | PERCENT
                           Only include codons containing a mutation present in
                           this percentage of sequences in the clone or more. 0
                           is all sequences. Converts the unincluded codons to
                           gaps.
  -l,--add-length          Whether to append the length of the sequence to the
                           end of the header, calculated after
                           --convert-to-amino-acids if enabled
  -N,--remove-N            Whether to replace N or n in the sequence with a gap,
                           '-'
  -g,--remove-germlines    Whether to remove germlines.
  -h,--remove-highly-mutated
                           Whether to remove highly mutated clone sequences (a
                           third of their sequence are different amino acids).
  -s,--remove-stops        Whether to remove sequences with stop codons
  -d,--legacy-remove-duplicates
                           Whether to remove duplicate sequences. LEGACY ONLY
  -O,--remove-out-of-frame Whether to remove sequences that are out of frame--if
                           the sequences or number of gaps is not divisible by 3
  -Y,--remove-unknown-nucleotides
                           Convert unknown nucleotides (not ACGTN-.) to gaps (-)
  --input-inframe [ ] | FIELD
                           Represents the 1 indexed field split by '|'
                           containing the inframe value (frames are 0, 1, or 2
                           like USCS definitions). For use with trim-frame.
  --input-outframe [ ] | FIELD
                           Represents the 1 indexed field split by '|'
                           containing the outframe value (frames are 0, 1, or 2
                           like USCS definitions). For use with trim-frame.
  -y,--trim-frame          Trim each sequence to be in frame by remove extra
                           nucleotides at the end. If input-inframe or
                           input-outframe is specified, follow those rules
                           instead.
  -r,--input-stop-range [106]|INT
                           Only search for stops with remove-stops up to this
                           amino acid position
  -c,--input-codon-mut [-1]|0|1|2|3
                           Only include codons with this many mutations or less
                           or more, depending on input-codon-mut-type (-1 is the
                           same as include all codons). Converts the unincluded
                           codon to gaps.
  -T,--input-codon-mut-type [=]|>|<
                           Only include codons with this many mutations (=) (or
                           lesser (<) or greater (>), depending on
                           input-codon-mut). Converts the unincluded codon to
                           gaps.
  -M,--input-mut-type [All]|Silent|Replacement
                           Only include codons with this all mutations (All),
                           (or silent (Silent) or replacement (Replacement)).
  -e,--input-change-field ((FIELD (Int), VALUE (String))
                           Change a field to a match, so a regex "ch.*_" to
                           field 2 of ">abc|brie_cheese_dude" would result in
                           ">abc|cheese_". Useful for getting specific
                           properties from a field. Can take a list of format
                           "(Int, String)&&(Int, String)&& ..." and so on. The
                           String is in regex format (POSIX extended). The first
                           in the tuple is the location of the field (1 indexed,
                           split by '|').
  -f,--input-custom-filter ((FIELD (Int), VALUE (String))
                           A custom filter. Can take a list of format "(Int,
                           String)&&(Int, String)&& ..." and so on. The String
                           is in regex format (POSIX extended), so if the entire
                           string 'a' is in the whole field, then you need to
                           input '^a$' for the beginning to the end! The first
                           in the tuple is the location of the field (1 indexed,
                           split by '|'). If you want to apply to the entire
                           header, either have the location as 0 or exclude the
                           location altogether (, Day 3|IGHV3) for instance will
                           match if the entire header is '>Day 3|IGHV3'. This
                           list will be filtered one at a time, so you cannot
                           get multiple filters, but you can remove multiple
                           filters.
  -G,--legacy-custom-germline
                           Whether to apply the custom filter to germlines (>>)
                           instead of sequences (>). LEGACY ONLY
  -m,--custom-remove       Whether to remove the sequences containing the custom
                           filter as opposed to remove the sequences that don't
                           contain the filter
  -V,--legacy-gene-allele-field [1]|INT
                           The field (1 indexed) of the gene allele name. LEGACY
                           ONLY
  -v,--legacy-count        Do not save output, just count genes and alleles from
                           the results. Requires gene-allele-field. LEGACY ONLY
  -o,--output FILE         The output fasta file
```
