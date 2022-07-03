# base

Library for storing and rudimentary manipulation of nucleic acid and protein sequences built with Python3.

It serves to be better than `biopython` as in it is aware of each monomer type.

Stores nucleotides as uint8 (allows for 256 monomer types):

- Monomer Byte Allocation
   - ATCG is the order of nucleotides

| Number  | Data Type     |
|---------|---------------|
| 000-049 | DNA's         |
| 050-099 | RNA's         |
| 100-199 | Amino Acids   |
| 200     | Invalid DNA   |
| 201     | Ambiguous DNA |
| 210     | Invalid RNA   |
| 211     | Ambiguous RNA |
| 220     | Invalid AA    |
| 221     | Ambiguous AA  |
| 256     | Null          |

- Models:
    - Monomer:
        - stores number
        - attributes:
            - static
            - molecular weight
                - lookup table
                - expand array dtype to float
    - BaseStrand:
        - children:
            - Strand
                - array of monomers
            - Segment
                - start & end indexes
        - attributes:
            - type
        - methods:
            - molecular weight:
                - replace all numbers in sequence
                - sum across array
            - length
            - GC-ratio
            - reverse complement
            - transcribe
            - translate
    - BaseSequence:
        - n-d array of strands
        - methods:
            - molecular weight
        - Children:
            - Protein Complex:
                - named array of n-dim of protein strands simulating protein complex
            - Sequence
                - 2-d array of monomers
                - children:
                    - Chromosome
                        - n-d array of monomers representing n-ploids
                    - CircularSequence
                        - methods:
                            - origin of replication
                            - shifting starting index