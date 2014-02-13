modify-fasta
============

**Gregory W. Schwartz**

```
Usage: modify_fasta [-i|--input FILE] [-a|--aminoAcids] [-A|--normalFasta] [-C|--convertToAminoAcids] [-N|--removeN] [-g|--removeGermlines] [-h|--removeHighlyMutated] [-s|--removeStops] [-r|--inputStopRange [106]|INT] [-c|--inputCodonMut [-1]|0|1|2|3] [-T|--inputCodonMutType [=]|>|<] [-M|--inputMutType [All]|Silent|Replacement] [-f|--inputCustomFilter ((FIELD_LOCATION (Int), FIELD_VALUE (String))] [-I|--infixCustomFilter] [-G|--customGermline] [-m|--customRemove] [-o|--output FILE]
  Modify fasta (and CLIP) files in several optional ways

Available options:
  -h,--help                Show this help text
  -i,--input FILE          The input CLIP fasta file
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
  -f,--inputCustomFilter ((FIELD_LOCATION (Int), FIELD_VALUE (String)) A custom filter. Can take a list of format "(Int, String)&&(Int, String)&& ..." and so on. The first in the tuple is the location of the field (1 indexed, split by '|'). If you want to apply to the entire header, either have the location as 0 or exclude the location altogether (, Day 3|IGHV3) for instance will match if the entire header is '>Day 3|IGHV3'. This list will be filtered one at a time, so you cannot get multiple filters, but you can remove multiple filters.
  -I,--infixCustomFilter   Whether to find the custom filter in the field as an infix (some part of the custom field matches the header) as opposed to an exact match (the entire field must be the custom field)
  -G,--customGermline      Whether to apply the custom filter to germlines (>>) instead of sequences (>)
  -m,--customRemove        Whether to remove the sequences containing the custom filter as opposed to remove the sequences that don't contain the filter
  -o,--output FILE         The output fasta file
```
