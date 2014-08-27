modify-fasta
============

**Gregory W. Schwartz**

This program will take a fasta file, do certain filterings and transformations
as specified by the user, and return the new fasta file.

To install:
```
cabal update
cabal install
```

```
Modify Fasta, Gregory W. Schwartz

Usage: modify-fasta [-i|--input FILE] [-a|--amino-acids] [-A|--clip-fasta]
                    [-C|--convert-to-amino-acids] [-N|--remove-N]
                    [-g|--remove-germlines] [-h|--remove-highly-mutated]
                    [-s|--remove-stops] [-d|--remove-duplicates]
                    [-r|--input-stop-range [106]|INT]
                    [-c|--input-codon-mut [-1]|0|1|2|3]
                    [-T|--input-codon-mut-type [=]|>|<]
                    [-M|--input-mut-type [All]|Silent|Replacement]
                    [-f|--input-custom-filter ((FIELD_LOCATION (Int), FIELD_VALUE (String))]
                    [-I|--infix-custom-filter] [-G|--custom-germline]
                    [-m|--custom-remove] [-V|--gene-allele-field [1]|INT]
                    [-v|--count] [-o|--output FILE]
  Modify fasta (and CLIP) files in several optional ways

Available options:
  -h,--help                Show this help text
  -i,--input FILE          The input fasta file or CLIP fasta file
  -a,--amino-acids         Whether these sequences are composed of amino acids
                           (AminoAcid) or nucleotides (Nucleotide)
  -A,--clip-fasta          Whether the input is a clip fasta file (has germline
                           >> sequences)
  -C,--convert-to-amino-acids
                           Whether to convert the filtered sequences to amino
                           acids in the output
  -N,--remove-N            Whether to remove N or n in the sequence
  -g,--remove-germlines    Whether to remove germlines
  -h,--remove-highly-mutated
                           Whether to remove highly mutated clone sequences (a
                           third of their sequence are different amino acids)
  -s,--remove-stops        Whether to remove sequences with stop codons
  -d,--remove-duplicates   Whether to remove duplicate sequences
  -r,--input-stop-range [106]|INT
                           Only search for stops with remove-stops up to this
                           amino acid position
  -c,--input-codon-mut [-1]|0|1|2|3
                           Only include codons with this many mutations or less
                           or more, depending on input-codon-mut-type (-1 is the
                           same as include all codons). Converts the codon to
                           gaps
  -T,--input-codon-mut-type [=]|>|<
                           Only include codons with this many mutations (=) (or
                           lesser (<) or greater (>), depending on
                           input-codon-mut). Converts the codon to gaps
  -M,--input-mut-type [All]|Silent|Replacement
                           Only include codons with this all mutations (All),
                           (or silent (Silent) or replacement (Replacement))
  -f,--input-custom-filter ((FIELD_LOCATION (Int), FIELD_VALUE (String))
                           A custom filter. Can take a list of format "(Int,
                           String)&&(Int, String)&& ..." and so on. The first in
                           the tuple is the location of the field (1 indexed,
                           split by '|'). If you want to apply to the entire
                           header, either have the location as 0 or exclude the
                           location altogether (, Day 3|IGHV3) for instance will
                           match if the entire header is '>Day 3|IGHV3'. This
                           list will be filtered one at a time, so you cannot
                           get multiple filters, but you can remove multiple
                           filters.
  -I,--infix-custom-filter Whether to find the custom filter in the field as an
                           infix (some part of the custom field matches the
                           header) as opposed to an exact match (the entire
                           field must be the custom field)
  -G,--custom-germline     Whether to apply the custom filter to germlines (>>)
                           instead of sequences (>)
  -m,--custom-remove       Whether to remove the sequences containing the custom
                           filter as opposed to remove the sequences that don't
                           contain the filter
  -V,--gene-allele-field [1]|INT
                           The field (1 indexed) of the gene allele name
  -v,--count               Do not save output, just count genes and alleles from
                           the results. Requires gene-allele-field
  -o,--output FILE         The output fasta file
```
