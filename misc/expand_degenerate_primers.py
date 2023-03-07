import sys

_degenerate_base_codes_ = {
    "R" : ["A","G"],
    "Y" : ["C","T"],
    "M" : ["A","C"],
    "K" : ["G","T"],
    "S" : ["G","C"],
    "W" : ["A","T"],
    "H" : ["A","T","C"],
    "B" : ["G","T","C"],
    "D" : ["G","A","T"],
    "V" : ["G","A","C"],
    "N" : ["A","C","G","T"]
}

def expand_primer(primer):
    ''' Expand degenerate bases in a primer
    :param str primer
    :yields primer(s) with degenerate bases expanded
    '''
    assert primer!="", "Require valid primer sequence!"

    degen_bases = []
    degen_bases_loc = []
    for base in _degenerate_base_codes_:
        loc = primer.find(base)
        if loc != -1:
            degen_bases.append(_degenerate_base_codes_[base])
            degen_bases_loc.append(loc)

    if len(degen_bases) == 0:
        yield primer
    else:
        combinations = [[]]
        for b in degen_bases:
            combinations = [i + [y] for y in b for i in combinations]

        for bases in combinations:
            temp = list(primer)
            for i,base in enumerate(bases):
                base_to_modify = degen_bases_loc[i]
                assert temp[base_to_modify] not in ["A","C","G","T"], "Error in Logic !!"
                temp[base_to_modify] = base

            yield "".join(temp)


def iterate_primer_file(primer_file):
    ''' Iterate over the primer file , expand degenerate bases and write to a new file
    :param str primer_file : A tsv file with information on the primers
    '''
    with open(primer_file,"r") as IN, open(primer_file+".expanded","w") as OUT:
        for line in IN:
            contents = line.strip("\n\r").split("\t")
            amplicon_id,forward_p_id,forward_p,reverse_p_id,reverse_p = contents

            if line.startswith("Region"): # header
                out = [amplicon_id,forward_p_id,"Forward Primer Name",forward_p,reverse_p_id,"Reverse Primer Name",reverse_p]
                OUT.write("\t".join(out)+"\n")
                continue

            if forward_p != "":
                forward_expanded_primers = list(expand_primer(forward_p.upper()))
            else:
                forward_expanded_primers = [""]

            if reverse_p != "":
                reverse_expanded_primers = list(expand_primer(reverse_p.upper()))
            else:
                reverse_expanded_primers = [""]

            i=1
            j=1
            for pf in forward_expanded_primers:
                for pr in reverse_expanded_primers:
                    out = [amplicon_id,forward_p_id,forward_p_id+"_"+str(i),pf,reverse_p_id,reverse_p_id+"_"+str(j),pr]
                    OUT.write("\t".join(out)+"\n")
                    j+=1
                i+=1

iterate_primer_file(sys.argv[1])
