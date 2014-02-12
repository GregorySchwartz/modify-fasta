modify-fasta
============

Modify Fasta, Gregory W. Schwartz

```
Usage: modify_fasta [-i|--input FILE] [-a|--aminoAcids] [-A|--normalFasta] [-C|--convertToAminoAcids] [-N|--removeN] [-g|--removeGermlines] [-h|--removeHighlyMutated] [-s|--removeStops] [-r|--inputStopRange [106]|INT] [-c|--inputCodonMut [-1]|0|1|2|3] [-T|--inputCodonMutType [=]|>|<] [-M|--inputMutType [All]|Silent|Replacement] [-f|--inputCustomFilter FIELD_VALUE (String)] [-F|--inputCustomField FIELD_LOCATION (Int)] [-I|--infixCustomFilter] [-G|--customGermline] [-m|--customRemove] [-o|--output FILE]

  Modify fasta (and CLIP) files in several optional ways

  Available options:

    -h,--help Show this help text
    -i,--input FILE      The input CLIP fasta file
    -a,--aminoAcids          Whether these sequences are composed of amino acids (AminoAcid) or nucleotides (Nucleotide)
    -A,--normalFasta         Whether the input is a normal fasta file (no germline >> sequences)
    -C,--convertToAminoAcids Whether to convert the filtered sequences to amino acids in the output
    -N,--removeN             Whether to remove N or n in the sequence
    -g,--removeGermlines     Whether to remove germlines
    -h,--removeHighlyMutated Whether to remove highly mutated clone sequences (a third of their sequence are different amino acids)
    -s,--removeStops         Whether to remove clone sequences with stop codons
    -r,--inputStopRange [106]|INT Only search for stops with removeStops up to this amino acid position
    -c,--inputCodonMut [-1]|0|1|2|3 Only include codons with this many mutations or less or more, depending on inputCodonMutType (-1 is the same as include all codons). Converts the codon to gaps
    -T,--inputCodonMutType [=]|>|< Only include codons with this many mutations (=) (or lesser (<) or greater (>), depending on inputCodonMut). Converts the codon to gaps
    -M,--inputMutType [All]|Silent|Replacement Only include codons with this all mutations (All), (or silent (Silent) or replacement (Replacement))
    -f,--inputCustomFilter FIELD_VALUE (String) A custom field value to keep, discarding all others
    -F,--inputCustomField FIELD_LOCATION (Int) The location of the field in inputCustomFilter (split by |, 0 means search the whole header)
    -I,--infixCustomFilter   Whether to find the custom filter in the field as an infix (some part of the custom field matches the header) as opposed to an exact match (the entire field must be the custom field)
    -G,--customGermline      Whether to apply the custom filter to germlines (>>) instead of sequences (>)
    -m,--customRemove        Whether to remove the sequences containing the custom filter as opposed to remove the sequences that don't contain the filter
    -o,--output FILE         The output fasta file
```
