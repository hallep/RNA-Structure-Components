# RNA-Structure-Components
Get and parse RNA secondary structure component information from bpRNA (https://github.com/hendrixlab/bpRNA)

### Install Requirements
* ct2db
* perl
* pandas
* numpy

## Code
All python functions can be found in ```RNAStructureComponents/components.py```

To run, call ```get_structure()```, passing in the file path to the .ct file of the target secondary structure

Alternatively, call ```run_bpRNA()```, passing in the file path to the .dbn file of the target secondary structure <br>
Then pass the output to ```parse_st()```

### Return Value
**```st_info```**: dictionary containing information about the different secondary structure components

#### Keys & Values
**"seq":** _str_ <br>
  the nucleotide sequence of the secondary structure

**"dbn":** _str_ <br>
  the dot-bracket notation of the secondary structure

**"sa":** _str_ <br>
  structure array: the position-wise component identity of each base

**"nk":** _str_ <br>
  knot annotation: a position-wise flag indicating non-canonical structure components <br>
  "N" indicates that the base is part of a canonical structure component <br>
  "K" indicates that the base is part of a pseudoknot or other non-canonical structure component 

**"ids":** _pandas.DataFrame_ <br>
  position-wise component identities <br>
  each column ("S", "B", "H", "I", "M", "PK", "X", "E", "NCBP", "segs", "id") is a different structure component <br>
  row indices are the 1-indexed positions in the structure <br>
  a value of -1 indicates that the base at position _{row}_ is not of component _{column}_ <br>
  a value of $\ge$ 0 indicates the index of the associated _{column}_ component for base _{row}_ in the corresponding pandas.DataFrame <br>

**"stems":** _pandas.DataFrame_ <br>
  stem component information <br>
  <ins>keys</ins>: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len" <br>
  see the ```parse_stems()``` subfunction within ```parse_st()``` for more information

**"bulges":** _pandas.DataFrame_ <br>
  bulge component information <br>
  <ins>keys</ins>: "rng", "coord", "seq", "len" <br>
  see the ```parse_bulges()``` subfunction within ```parse_st()``` for more information

**"hairpins":** _pandas.DataFrame_ <br>
  hairpin loop component information <br>
  <ins>keys</ins>: "rng", "coord", "seq", "len" <br>
  see the ```parse_hairpins()``` subfunction within ```parse_st()``` for more information

**"iloops":** _pandas.DataFrame_ <br>
  internal loop component information <br>
  <ins>keys</ins>: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len5", "len3" <br>
  see the ```parse_iloops()``` subfunction within ```parse_st()``` for more information

**"mloops":** _pandas.DataFrame_ <br>
  multiloop component information <br>
  <ins>keys</ins>: "rng", "coord", "seq", "len" <br>
  see the ```parse_multiloops()``` subfunction within ```parse_st()``` for more information

**"pknots":** _pandas.DataFrame_ <br>
  pseudoknot component information <br>
  <ins>keys</ins>: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len" <br>
  see the ```parse_pseudoknots()``` subfunction within ```parse_st()``` for more information

**"xloops":** _pandas.DataFrame_ <br>
  external loop component information <br>
  <ins>keys</ins>: "rng", "coord", "seq", "len" <br>
  see the ```parse_xloops()``` subfunction within ```parse_st()``` for more information

**"ends":** _pandas.DataFrame_ <br>
  dangling end component information <br>
  <ins>keys</ins>: "rng", "coord", "seq", "len" <br>
  see the ```parse_ends()``` subfunction within ```parse_st()``` for more information

**"ncbp":** _pandas.DataFrame_ <br>
  non-canonical base pairing information <br>
  <ins>keys</ins>: "coord5", "coord3", "seq5", "seq3" <br>
  see the ```parse_ncbp()``` subfunction within ```parse_st()``` for more information

**"segs":** _pandas.DataFrame_ <br>
  segment information <br>
  <ins>keys</ins>: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len5", "len3" <br>
  see the ```parse_segments()``` subfunction within ```parse_st()``` for more information
