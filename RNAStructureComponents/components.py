import os
import subprocess
import pandas as pd
import numpy as np

struct_dict = {
    "S" : "stems",
    "H" : "hairpins",
    "B" : "bulges",
    "I" : "iloops",
    "M" : "mloops",
    "PK" : "pknots",
    "X" : "xloops",
    "E" : "ends"
}

def create_dbn(src:str, dst:str="bpRNA/structure.dbn"):
    
    ''' Convert .ct file into a .dbn that bpRNA can parse

    Parameters
    ----------
    src : str
        filepath to .ct file (relative to current directory)
    dst : str (default = bpRNA/structure.dbn)
        filepath to new .dbn file (relative to current directory)
    '''
    
    # create .dbn
    out = subprocess.check_output(f"ct2db {src}", shell=True, text=True)
    lines = out.split("\n")

    # define components
    header = ">" + src.split("/")[-1].split(".")[0]
    seq = lines[1]
    struct = lines[2]

    # write new .dbn file
    with open(dst, "w") as dbn:
        dbn.write(f"{header}\n{seq}\n{struct}")
    
def run_bpRNA(src:str="structure.dbn") -> str:

    ''' Run bpRNA

    Parameters
    ----------
    src : str (default = "structure.dbn")
        filepath to .dbn (or .bpseq) file (relative to {bpRNA} sub-directory) to pass to bpRNA

    Returns
    -------
    _ : str
        filepath to output .st file
    '''

    # call bpRNA perl script
    os.chdir("bpRNA")
    os.system(f"perl bpRNA.pl {src}")
    os.chdir("..")

    # get filepath to .st output file
    return f"bpRNA/{src.split("/")[-1].split(".")[0]}.st"

def parse_st(src:str="bpRNA/structure.st") -> dict[str, any]:

    ''' Parse .st file output from bpRNA
    
    Parameters
    ----------
    src : str (default = "bpRNA/structure.st")
        filepath to .st file (relative to current directory)
    
    Returns
    -------
    st_info : dict
        structure information
        * "seq" : str \\
            base sequence
        * "dbn" : str \\
            dot-bracket notation of structure
        * "sa" : str \\
            structure array - position-wise component identity of each base
        * "nk" : str \\
            knot annotation - position-wise flag indicating non-canonical structures \\
            "N": part of a canonical structure \\
            "K": part of a pseudoknot or other non-canonical structure \\
        * "ids" : pandas.DataFrame \\
            position-wise component identities \\
            columns : "S", "B", "H", "I", "M", "PK", "X", "E", "NCBP", "id" \\
            rows: 1-index position in structure \\
            (see below for details)
        * "stems" : pandas.DataFrame \\
            stem components information from parse_stems() \\
            keys: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len"
        * "bulges" : pandas.DataFrame \\
            bulge component information from parse_bulges() \\
            keys: "rng", "coord", "seq", "len"
        * "hairpins" : pandas.DataFrame \\
            hairpin loop component information from parse_hairpins() \\
            keys: "rng", "coord", "seq", "len"
        * "iloops" : pandas.DataFrame \\
            internal loop component information from parse_iloops() \\
            keys: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len5", "len3"
        * "mloops" : pandas.DataFrame \\
            multiloop component information from parse_multiloops() \\
            keys: "rng", "coord", "seq", "lengths", "len"
        * "pknots" : pandas.DataFrame \\
            pseudoknot component information from parse_pseudoknots() \\
            keys: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len"
        * "xloops" : pandas.DataFrame \\
            external loop component information from parse_xloops() \\
            keys: "rng", "coord", "seq", "len"
        * "ends" : pandas.DataFrame \\
            dangling end component information from parse_ends() \\
            keys: "rng", "coord", "seq", "len"
        * "ncbp" : pandas.DataFrame \\
            non-canonical base pairing information from parse_ncbp() \\
            keys: "coord5", "coord3", "seq5", "seq3"
        * "segs" : pandas.DataFrame \\
            segment information from parse_segments() \\
            keys: "rng5", "rng3", "coord5", "coord3", "seq5", "seq3", "len5", "len3"
    
    "ids" DataFrame
    ---------------
     * "S": int
        * -1: not part of stem
        * \\>=0: index in st_info["stems"]
     * "B": int
        * -1: not part of bulge
        * \\>=0: index in st_info["bulges"]
     * "H": int
        * -1: not part of hairpin loop
        * \\>=0: index in st_info["hairpins"]
     * "I": int
        * -1: not part of internal loop
        * \\>=0: index in st_info["iloops"]
     * "M": int
        * -1: not part of multiloop
        * \\>=0: index in st_info["mloops"]
     * "PK": int
        * -1: not part of pseudoknot
        * \\>=0: index in st_info["pknots"]
     * "X": int
        * -1: not part of external loop
        * \\>=0: index in st_info["xloops"]
     * "E": int
        * -1: not part of danging end
        * \\>=0: index in st_info["ends"]
     * "NCBP": int
        * -1: unpaired or canonically base paired
        * \\>=0: index in st_info["ncbp"]
     * "segs" : int
        * -1: not part of a segment
        * \\>=0: index in st_info["segs"]
     * "id": list of strings \\
        component ("S", "H", "B", "I", "M", "PK", "X", "E") pertaining to each base
    '''

    # read .st file
    with open(src, "r") as st:
        lines = st.readlines()

    header = [l for l in lines if l.startswith("#")]
    sequence = lines[len(header)].strip()
    dbn = lines[len(header)+1].strip()
    
    # sa = structure array - position-wise key identity of each base (see struct_dict)
    sa = lines[len(header)+2].strip()

    # nk = whether bases are canonically paired/unpaired (N) or non-canonically paired (K)
    nk = lines[len(header)+3].strip()

    # structural components
    components = lines[len(header)+4:]

    # position-wise component identities
    comp_ids = {
        "S" : [-1] * len(sequence),
        "B" : [-1] * len(sequence),
        "H" : [-1] * len(sequence),
        "I" : [-1] * len(sequence),
        "M" : [-1] * len(sequence),
        "PK" : [-1] * len(sequence),
        "X" : [-1] * len(sequence),
        "E" : [-1] * len(sequence),
        "NCBP" : [-1] * len(sequence),
        "segs" : [-1] * len(sequence),
        "id" : [[] for _ in range(len(sequence))]
    }

    # STEMS
    def parse_stems(stems:list) -> pd.DataFrame:

        ''' Parse stem components
        
        Parameters
        ----------
        stems : list
            list of strings containing stem information \\
            format: S# #..# "NNNNN" #..# "NNNNN"
        
        Returns
        -------
        df : pandas.DataFrame
            stem information \\
            note: coord1[i] pairs with coord2[i] and seq1[i] pairs with seq2[i]
            * "rng5" : tuple \\
                rng5[0] = 5' 1-indexed position of 5' side of stem \\
                rng5[1] = 3' 1-indexed position of 5' side of stem
            * "rng3" : tuple \\
                rng3[0] = 5' 1-indexed position of 3' side of stem \\
                rng3[1] = 3' 1-indexed position of 3' side of stem
            * "coord5" : numpy.NDArray \\
                1-indexed positions of 5' side of stem
            * "coord3" : numpy.NDArray \\
                1-indexed positions of 3' side of stem (in reverse order)
            * "seq5" : str \\
                bases on 5' side of stem
            * "seq3" : str \\
                bases on 3' side of stem (in reverse order)
            * "len" : int \\
                number of base pairs in stem
        '''

        ind = [0] * len(stems)
        rng5 = [0] * len(stems)
        rng3 = [0] * len(stems)
        coord5 = [0] * len(stems)
        coord3 = [0] * len(stems)
        seq5 = [""] * len(stems)
        seq3 = [""] * len(stems)
        length = [0] * len(stems)

        # for each stem:
        for i,s in enumerate(stems):

            # stem number
            ind[i] = int(s[0][1:])

            # 5' side coordinates
            r5 = [int(c) for c in s[1].split("..")]
            rng5[i] = tuple(r5)
            coord5[i] = np.arange(r5[0], r5[1]+1)

            # mark in comp_id dictionary
            for z in coord5[i]:
                comp_ids["S"][z-1] = ind[i]
                comp_ids["id"][z-1].append("S")
            
            # 5' side base sequence
            seq5[i] = s[2].strip("\"")

            # 3' side coordinates
            r3 = [int(c) for c in s[3].split("..")]
            rng3[i] = tuple(r3)
            coord3[i] = np.arange(r3[0], r3[1]+1)[::-1]
            
            # mark in comp_id dictionary
            for z in coord3[i]:
                comp_ids["S"][z-1] = ind[i]

                if "S" not in comp_ids["id"][z-1]:
                    comp_ids["id"][z-1].append("S")

            # 3' side base sequence
            seq3[i] = s[4].strip("\"")[::-1]

            # length
            length[i] = r5[1] - r5[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng5" : rng5,
            "rng3" : rng3,
            "coord5" : coord5,
            "coord3" : coord3,
            "seq5" : seq5,
            "seq3" : seq3,
            "len" : length
        }, index=ind)

        return df

    stems = parse_stems([c.split() for c in components if c.startswith("S")])

    # BULGES
    def parse_bulges(bulges:list) -> pd.DataFrame:

        ''' Parse bulge components
        
        Parameters
        ----------
        bulges : list
            list of strings containing bulge information \\
            format: B# #..# "NNNNN" (#,#) N:N (#,#) N:N
        
        Returns
        -------
        df : pandas.DataFrame
            bulge information
            * "rng" : tuple \\
                rng[0] = 5' 1-indexed position of bulge \\
                rng[1] = 3' 1-indexed position of bulge
            * "coord" : numpy.NDArray \\
                1-indexed positions of bulge
            * "seq" : str \\
                bases in bulge
            * "len" : int \\
                number of bases in bulge
        '''

        ind = [0] * len(bulges)
        rng = [0] * len(bulges)
        coord = [0] * len(bulges)
        seq = [""] * len(bulges)
        length = [0] * len(bulges)

        # for each bulge:
        for i,b in enumerate(bulges):

            # bulge number
            ind[i] = int(b[0][1:])

            # coordinates
            r = [int(c) for c in b[1].split("..")]
            rng[i] = tuple(r)
            coord[i] = np.arange(r[0], r[1]+1)

            # mark in comp_id dictionary
            for z in coord[i]:
                comp_ids["B"][z-1] = ind[i]
                comp_ids["id"][z-1].append("B")

            # base sequence
            seq[i] = b[2].strip("\"")

            # length
            length[i] = r[1] - r[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng" : rng,
            "coord" : coord,
            "seq" : seq,
            "len" : length
        }, index=ind)

        return df

    bulges = parse_bulges([c.split() for c in components if c.startswith("B")])

    # HAIRPIN LOOPS
    def parse_hairpins(hairpins:list) -> pd.DataFrame:

        ''' Parse hairpin loop components
        
        Parameters
        ----------
        hairpins : list
            list of strings containing hairpin loop information \\
            format: H# #..# "NNNNN" (#,#) N:N
        
        Returns
        -------
        df : pandas.DataFrame
            stem information
            * "rng" : tuple \\
                rng[0] = 5' 1-indexed position of hairpin loop \\
                rng[1] = 3' 1-indexed position of hairpin loop
            * "coord" : numpy.NDArray \\
                1-indexed positions of hairpin loop
            * "seq" : str \\
                bases in hairpin loop
            * "len" : int \\
                number of bases in hairpin loop
        '''

        ind = [0] * len(hairpins)
        rng = [0] * len(hairpins)
        coord = [0] * len(hairpins)
        seq = [""] * len(hairpins)
        length = [0] * len(hairpins)

        # for each hairpin loop:
        for i,h in enumerate(hairpins):

            # hairpin number
            ind[i] = int(h[0][1:])

            # coordinates
            r = [int(c) for c in h[1].split("..")]
            r.sort()
            rng[i] = tuple(r)
            coord[i] = np.arange(r[0], r[1]+1)
            
            # mark in comp_id dictionary
            for z in coord[i]:
                comp_ids["H"][z-1] = ind[i]
                comp_ids["id"][z-1].append("H")

            # base sequence
            if h[2] == "\"\"":
                seq[i] = sequence[r[0]-1:r[1]]
            else:
                seq[i] = h[2].strip("\"")

            # length
            length[i] = r[1] - r[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng" : rng,
            "coord" : coord,
            "seq" : seq,
            "len" : length
        }, index=ind)

        return df

    hairpins = parse_hairpins([c.split() for c in components if c.startswith("H")])

    # INTERNAL LOOPS
    def parse_iloops(iloops:list) -> pd.DataFrame:

        ''' Parse internal loop components
        
        Parameters
        ----------
        iloops : list
            list of strings containing internal loop information \\
            format: I#.# #..# "NNNNN" (#,#) N:N
        
        Returns
        -------
        df : pandas.DataFrame
            internal loop information
            * "rng5" : tuple \\
                rng5[0] = 5' 1-indexed position of 5' side of internal loop \\
                rng5[1] = 3' 1-indexed position of 5' side of internal loop
            * "rng3" : tuple \\
                rng3[0] = 5' 1-indexed position of 3' side of internal loop \\
                rng3[1] = 3' 1-indexed position of 3' side of internal loop
            * "coord5" : numpy.NDArray \\
                1-indexed positions of 5' side of internal loop
            * "coord3" : numpy.NDArray \\
                1-indexed positions of 3' side of internal loop
            * "seq5" : str \\
                bases on 5' side of internal loop
            * "seq3" : str \\
                bases on 3' side of internal loop
            * "len5" : int \\
                number of bases on 5' side of internal loop
            * "len3" : int \\
                number of bases on 3' side of internal loop
        '''

        ind = [0] * int(len(iloops) / 2)
        rng = [[0] * int(len(iloops) / 2) for _ in range(2)]
        coord = [[0] * int(len(iloops) / 2) for _ in range(2)]
        seq = [[""] * int(len(iloops) / 2) for _ in range(2)]
        length = [[0] * int(len(iloops) / 2) for _ in range(2)]

        # for each internal loop:
        for l in iloops:
            loop = int(l[0][1:].split(".")[0])
            i = loop - 1

            # side (5' or 3')
            s = 0 if l[0][1:].split(".")[1] == "1" else 1
            
            # internal loop number
            ind[i] = loop

            # coordinates
            r = [int(c) for c in l[1].split("..")]
            rng[s][i] = tuple(r)
            coord[s][i] = np.arange(r[0], r[1]+1)
            
            # mark in comp_id dictionary
            for z in coord[s][i]:
                comp_ids["I"][z-1] = loop
                comp_ids["id"][z-1].append("I")

            # base sequence
            seq[s][i] = l[2].strip("\"")

            # length
            length[s][i] = r[1] - r[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng5" : rng[0],
            "rng3" : rng[1],
            "coord5" : coord[0],
            "coord3" : coord[1],
            "seq5" : seq[0],
            "seq3" : seq[1],
            "len5" : length[0],
            "len3" : length[1]
        }, index=ind)

        return df

    iloops = parse_iloops([c.split() for c in components if c.startswith("I")])

    # MULTILOOPS
    def parse_multiloops(mloops:list) -> pd.DataFrame:

        ''' Parse multiloop components
        
        Parameters
        ----------
        multi : list
            list of strings containing multiloop information \\
            format: M# #..# "NNNNN" (#,#) N:N (#,#) N:N
        
        Returns
        -------
        df : pandas.DataFrame
            multiloop information
            * "rng" : list of tuples \\
                rng[i][0] = 5' 1-indexed position of multiloop sections \\
                rng[i][1] = 3' 1-indexed position of multiloop sections
            * "coord" : list of numpy.NDArray \\
                1-indexed positions of multiloop sections
            * "seq" : list of strings \\
                bases in multiloop sections
            * "lengths" : list of ints \\
                number of bases in multiloop sections
            * "len" : int \\
                total number of bases in multiloop
        '''
        
        multi = {}

        # for each multiloop:
        for m in mloops:

            # multiloop number
            i = int(m[0][1:].split(".")[0])

            if i not in multi.keys():
                multi[i] = [[], [], [], [], 0]

            # coordinates
            r = [int(c) for c in m[1].split("..")]
            r.sort()
            multi[i][0].append(tuple(r))
            multi[i][1].append(np.arange(r[0], r[1]+1))
            
            # mark in comp_id dictionary
            for z in multi[i][1][-1]:
                comp_ids["M"][z-1] = i
                comp_ids["id"][z-1].append("M")

            # base sequence
            if m[2] == "\"\"":
                multi[i][2].append(sequence[r[0]-1:r[1]])
            else:
                multi[i][2].append(m[2].strip("\""))

            # length
            l = r[1] - r[0] + 1
            multi[i][3].append(l)
            multi[i][4] += l

        # create DataFrame
        df = pd.DataFrame({
            "rng" : [m[0] for m in multi.values()],
            "coord" : [m[1] for m in multi.values()],
            "seq" : [m[2] for m in multi.values()],
            "lengths" : [m[3] for m in multi.values()],
            "len" : [m[4] for m in multi.values()]
        }, index=multi.keys())

        return df

    mloops = parse_multiloops([c.split() for c in components if c.startswith("M")])

    # PSEUDOKNOTS
    def parse_pseudoknots(pseudo:list) -> pd.DataFrame:

        ''' Parse pseudoknot components
        
        Parameters
        ----------
        pseudo : list
            list of strings representing pseudoknots \\
            format1: PK# #bp #..# #..# Z#.# #..# Z#.# #..# \\
            format2: PK#.# # N # N
        
        Returns
        -------
        df : pandas.DataFrame
            pseudoknot information \\
            note: coord1[i] pairs with coord2[i] and seq1[i] pairs with seq2[i]
            * "rng5" : tuple \\
                rng5[0] = 5' 1-indexed position of 5' side of pseudoknot \\
                rng5[1] = 3' 1-indexed position of 5' side of pseudoknot
            * "rng3" : tuple \\
                rng3[0] = 5' 1-indexed position of 3' side of pseudoknot \\
                rng3[1] = 3' 1-indexed position of 3' side of pseudoknot
            * "coord5" : numpy.NDArray \\
                1-indexed positions of 5' side of pseudoknot
            * "coord3" : numpy.NDArray \\
                1-indexed positions of 3' side of pseudoknot
            * "seq5" : str \\
                bases on 5' side of pseudoknot
            * "seq3" : str \\
                bases on 3' side of pseudoknot
            * "len" : int \\
                number of base pairs in pseudoknot
        '''
        
        pknots = {}

        # for each pseudoknot:
        for k in pseudo:
            
            # initial definition
            if "." not in k[0]:

                # pseudoknot number
                i = int(k[0][2:])

                # size
                s = int(k[1][:-2])

                # 5' range
                r5 = tuple([int(c) for c in k[2].split("..")])

                # 3' range
                r3 = tuple([int(c) for c in k[3].split("..")])

                pknots[i] = [r5] + [[0] * s for _ in range(2)] + [r3] + [[0] * s for _ in range(2)] + [s]
            
            # specifics
            else:

                # pseudoknot number
                i = int(k[0][2:].split(".")[0])

                # position number
                j = int(k[0][1:].split(".")[1]) - 1

                # 5' coordinates
                pknots[i][1][j] = int(k[1])

                # mark in comp_id dictionary
                comp_ids["PK"][int(k[1])-1] = i
                comp_ids["id"][int(k[1])-1].append("PK")
                
                # 5' base sequence
                pknots[i][2][j] = k[2]

                # 3' coordinates
                pknots[i][4][j] = int(k[3])
                
                # mark in comp_id dictionary
                comp_ids["PK"][int(k[3])-1] = i
                comp_ids["id"][int(k[3])-1].append("PK")

                # 3' base sequence
                pknots[i][5][j] = k[4]
        
        # create DataFrame
        df = pd.DataFrame({
            "rng5" : [k[0] for k in pknots.values()],
            "rng3" : [k[3] for k in pknots.values()],
            "coord5" : [k[1] for k in pknots.values()],
            "coord3" : [k[4] for k in pknots.values()],
            "seq5" : [k[2] for k in pknots.values()],
            "seq3" : [k[5] for k in pknots.values()],
            "len" : [k[6] for k in pknots.values()]
        }, index=pknots.keys())

        return df

    pknots = parse_pseudoknots([c.split() for c in components if c.startswith("PK")])

    # EXTERNAL LOOPS
    def parse_xloops(xloops:list) -> pd.DataFrame:

        ''' Parse external loop components
        
        Parameters
        ----------
        xloops : list
            list of strings containing external loop information \\
            format: X# #..# "NNNNN" (#,#) N:N (#,#) N:N
        
        Returns
        -------
        df : pandas.DataFrame
            external loop information
            * "rng" : tuple \\
                rng[0] = 5' 1-indexed position of external loop \\
                rng[1] = 3' 1-indexed position of external loop
            * "coord" : numpy.NDArray \\
                1-indexed positions of external loop
            * "seq" : str \\
                bases in external loop
            * "len" : int \\
                number of bases in external loop
        '''

        ind = [0] * len(xloops)
        rng = [0] * len(xloops)
        coord = [0] * len(xloops)
        seq = [""] * len(xloops)
        length = [0] * len(xloops)

        # for each external loop:
        for i,x in enumerate(xloops):

            # external loop number
            ind[i] = int(x[0][1:])

            # coordinates
            r = [int(c) for c in x[1].split("..")]
            rng[i] = tuple(r)
            coord[i] = np.arange(r[0], r[1]+1)

            # mark in comp_id dictionary
            for z in coord[i]:
                comp_ids["X"][z-1] = ind[i]
                comp_ids["id"][z-1].append("X")

            # base sequence
            seq[i] = x[2].strip("\"")

            # length
            length[i] = r[1] - r[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng" : rng,
            "coord" : coord,
            "seq" : seq,
            "len" : length
        }, index=ind)

        return df

    xloops = parse_xloops([c.split() for c in components if c.startswith("X")])

    # DANGLING ENDS
    def parse_ends(ends:list) -> pd.DataFrame:

        ''' Parse dangling end components
        
        Parameters
        ----------
        ends : list
            list of strings containingdangling end information \\
            format: E# #..# "NNNNN"
        
        Returns
        -------
        df : pandas.DataFrame
            dangling end information
            * "rng" : tuple \\
                rng[0] = 5' 1-indexed position of dangling end \\
                rng[1] = 3' 1-indexed position of dangling end
            * "coord" : numpy.NDArray \\
                1-indexed positions of danging end
            * "seq" : str \\
                bases in danging end
            * "len" : int \\
                number of bases in dangling end
        '''

        ind = [0] * len(ends)
        rng = [0] * len(ends)
        coord = [0] * len(ends)
        seq = [""] * len(ends)
        length = [0] * len(ends)

        # for each danging end:
        for i,e in enumerate(ends):

            # danging end number
            ind[i] = int(e[0][1:])

            # coordinates
            r = [int(c) for c in e[1].split("..")]
            rng[i] = tuple(r)
            coord[i] = np.arange(r[0], r[1]+1)

            # mark in comp_id dictionary
            for z in coord[i]:
                comp_ids["E"][z-1] = ind[i]
                comp_ids["id"][z-1].append("E")

            # base sequence
            seq[i] = e[2].strip("\"")

            # length
            length[i] = r[1] - r[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng" : rng,
            "coord" : coord,
            "seq" : seq,
            "len" : length
        }, index=ind)

        return df

    ends = parse_ends([c.split() for c in components if c.startswith("E")])

    # NON-CANONICAL BASE PAIRS
    def parse_ncbp(pairs:list) -> pd.DataFrame:

        ''' Parse non-canonical base pair components
        
        Parameters
        ----------
        pairs : list
            list of strings representing non-canonical base pairs \\
            format2: NCBP# # N # N
        
        Returns
        -------
        df : pandas.DataFrame
            non-canonical base pair information \\
            note: coord1[i] pairs with coord2[i] and seq1[i] pairs with seq2[i]
            * "coord5" : int \\
                1-indexed position of 5' side of pair
            * "coord3" : int \\
                1-indexed position of 3' side of pair
            * "seq5" : str \\
                base on 5' side of pair
            * "seq3" : str \\
                base on 3' side of pair
        '''
        
        ind = [0] * len(pairs)
        coord5 = [0] * len(pairs)
        coord3 = [0] * len(pairs)
        seq5 = [""] * len(pairs)
        seq3 = [""] * len(pairs)

        # for each base pair:
        for i,p in enumerate(pairs):
            
            # pair number
            ind[i] = int(p[0][4:])

            # 5' coordinates
            coord5[i] = int(p[1])

            # mark in comp_id dictionary
            comp_ids["NCBP"][int(p[1])-1] = ind[i]
            
            # 5' base sequence
            seq5[i] = p[2]

            # 3' coordinates
            coord3[i] = int(p[3])
            
            # mark in comp_id dictionary
            comp_ids["NCBP"][int(p[3])-1] = ind[i]

            # 3' base sequence
            seq3[i] = p[4]
        
        # create DataFrame
        df = pd.DataFrame({
            "coord5" : coord5,
            "coord3" : coord3,
            "seq5" : seq5,
            "seq3" : seq3
        }, index=ind)

        return df

    ncbps = parse_ncbp([c.split() for c in components if c.startswith("NCBP")])

    # SEGMENTS
    def parse_segments(segs:list) -> pd.DataFrame:

        ''' Parse segments
        
        Parameters
        ----------
        segs : list
            list of strings containing segment information \\
            format: segment# #bp #..# NNNNN #..# NNNNN
        
        Returns
        -------
        df : pandas.DataFrame
            segment information
            * "rng5" : tuple \\
                rng5[0] = 5' 1-indexed position of 5' side of segment \\
                rng5[1] = 3' 1-indexed position of 5' side of segment
            * "rng3" : tuple \\
                rng3[0] = 5' 1-indexed position of 3' side of segment \\
                rng3[1] = 3' 1-indexed position of 3' side of segment
            * "coord5" : numpy.NDArray \\
                1-indexed positions of 5' side of segment
            * "coord3" : numpy.NDArray \\
                1-indexed positions of 3' side of segment
            * "seq5" : str \\
                bases on 5' side of segment
            * "seq3" : str \\
                bases on 3' side of segment
            * "len5" : int \\
                number of bases on 5' side of segment
            * "len3" : int \\
                number of bases on 3' side of segment
        '''

        ind = [0] * len(segs)
        rng5 = [0] * len(segs)
        rng3 = [0] * len(segs)
        coord5 = [0] * len(segs)
        coord3 = [0] * len(segs)
        seq5 = [""] * len(segs)
        seq3 = [""] * len(segs)
        len5 = [0] * len(segs)
        len3 = [0] * len(segs)

        # for each stem:
        for i,s in enumerate(segs):

            # stem number
            ind[i] = int(s[0][7:])

            # 5' side coordinates
            r5 = [int(c) for c in s[2].split("..")]
            rng5[i] = tuple(r5)
            coord5[i] = np.arange(r5[0], r5[1]+1)

            # mark in comp_id dictionary
            for z in coord5[i]:
                comp_ids["segs"][z-1] = ind[i]

            # 5' side base sequence
            seq5[i] = s[3]

            # 5' length
            len5[i] = r5[1] - r5[0] + 1

            # 3' side coordinates
            r3 = [int(c) for c in s[4].split("..")]
            rng3[i] = tuple(r3)
            coord3[i] = np.arange(r3[0], r3[1]+1)[::-1]
            
            # mark in comp_id dictionary
            for z in coord3[i]:
                comp_ids["segs"][z-1] = ind[i]

            # 3' side base sequence
            seq3[i] = s[5]
            
            # 3' length
            len3[i] = r3[1] - r3[0] + 1

        # create DataFrame
        df = pd.DataFrame({
            "rng5" : rng5,
            "rng3" : rng3,
            "coord5" : coord5,
            "coord3" : coord3,
            "seq5" : seq5,
            "seq3" : seq3,
            "len5" : len5,
            "len3" : len3
        }, index=ind)

        return df

    segs = parse_segments([c.split() for c in components if c.startswith("segment")])

    # POSITION-WISE COMPONENT IDENTITIES
    ids = pd.DataFrame(comp_ids, index=[i+1 for i in range(len(sequence))])

    # create and return final dictionary
    st_info = {
        "seq" : sequence,
        "dbn" : dbn,
        "sa" : sa,
        "nk" : nk,
        "ids" : ids,
        "stems" : stems,
        "bulges" : bulges,
        "hairpins" : hairpins,
        "iloops" : iloops,
        "mloops" : mloops,
        "pknots" : pknots,
        "xloops" : xloops,
        "ends" : ends,
        "ncbp" : ncbps,
        "segs" : segs
    }

    return st_info

def get_structure(src:str) -> dict[str, any]:
    
    ''' Get RNA secondary structure information
    
    Parameters
    ----------
    src : str (default = None)
        filepath to .ct file containing secondary structure \\
        Note: {src} should only contain one .ct structure

    Returns
    -------
    st_info : dict
        dictionary of structure information from parse_st()
        * keys: "seq", "dbn", "nk", "ids", "stems", "bulges", "hairpins",
            "iloops", "mloops", "pknots", "xloops", "ends", "ncbp", "segs"
    '''

    # create dot-bracket notation file
    create_dbn(src)

    # run bpRNA
    st = run_bpRNA()

    # parse .st file
    st_info = parse_st(st)

    # return component information
    return st_info
