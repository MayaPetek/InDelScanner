def classify_dna(dna_error):
    """
    Given a dna_error tuple, it describes the type of mutation: 
    d(eletion), i(nsertion), s(ubstitution), f(rameshift), b(roken), wt 
    or a combination of these.
    """
    if dna_error is None:  # empty or broken reads
        return 'b'
    elif dna_error == ():
        return 'wt'
    elif len(dna_error) > 1:
        # expect substitutions
        if dna_error[-1][1] == 'f':  # frameshifts are always the last mutation
            return 'f'
        else:
            for k in range(len(dna_error)):
                if dna_error[k][1] == 'i':
                    return 'si'
                elif dna_error[k][1] == 'd':
                    return 'sd'
            return 's'
    elif len(dna_error) == 1:
        if dna_error[0][1] == 'f':
            return 'f'
        elif dna_error[0][1] == 'd':
            return 'd' + dna_error[0][2]  # length of deletion
        elif dna_error[0][1] == 'i':
            return 'i' + str(len(dna_error[0][2]))
        elif dna_error[0][1] == 's':
            return 's'
    else:
        return 'b'

def classify_protein(mutation):
    """
    Detect whether a protein mutation is an insertion/deletion/substituion and if more than 1, (non)consecutive
    :param mutation:
    :return:
    """
    if mutation is None:  # came from empty or broken reads
        return 'b'
    elif mutation == ():
        return 'b'
    else:
        m = []
        for pos in range(len(mutation)):
            t = mutation[pos][1]
            if t != 'i':
                m.append(t)
            else:
                m.append(t + str(len(mutation[pos][2])))
        # need to distinguish between consecutive mutations and likely sequencing errors
        if 'f' in m:
            return 'f'
        elif len(m) == 1:
            if mutation[0][-1] == '*':
                return '*'
            else:
                return ''.join(m)
        elif len(m) >= 1:
            if mutation[-1][-1] == '*':
                return '*'
            else:
                c = is_mutation_consecutive(mutation) + '-' + ''.join(m)
                return c


def is_mutation_consecutive(mutation):
    """
    If only consecutve amino acids are affected, return 'c', else return 'nc'
    :param mutation:
    :return:
    """
    for pos in range(1, len(mutation)):
        if mutation[pos][0] != (mutation[pos - 1][0] + 1) :
            return 'nc'
    return 'c'


def translate_variant_type(dna_type):
    """
    Given a type of dna mutation (d3, d6 etc.), it translates this
    into what effect this should have on protein for target mutations.
    """
    dna_to_prot_type = {'d3':['d', 'c-sd', '*'], 'd6':['c-dd', 'c-sdd', '*'], 'd9':['c-ddd', 'c-sddd', '*'],
                        'i3':['i1','c-si1', '*'], 'i6': ['i2','c-si2', '*'], 'i9':['i3','c-si3', '*'],
                        's': ['s', 'c-ss', '*']
                       }
    
    try:
        return dna_to_prot_type[dna_type]
    except KeyError:
        return []      

def get_dna_composition(all_references, cutoff=10):
    """
    Classify all detected mutations by what kind they are: substitution, deletion, insertion (including length), frameshift, or other.
    Collect the results in a dataframe, both for mutation counts and associated sequencing reads.
    Relatively simple logic because two different DNA mutations can give rise to one protein mutation, but not vice versa.
    """

    dna_count = {}
    dna_reads = {}
    
    dna_to_stop = {'d3': 0, 'd6': 0, 'd9': 0}
    
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            # print('Analysing background {0} and fraction {1}. Unusal mutations: '.format(background, fraction))
            distinct_mutations = 0
            total_count = 0
            dna_count[background + '.' + fraction] = {'s': 0, 'd3': 0, 'd6': 0, 'd9': 0, 'i3': 0, 'i6': 0, 'i9': 0, 'f': 0,
                                                'b': 0, 'other': 0, 'sd': 0, 'si': 0}
            dna_reads[background + '.' + fraction] = {'s': 0, 'd3': 0, 'd6': 0, 'd9': 0, 'i3': 0, 'i6': 0, 'i9': 0, 'f': 0,
                                                'b': 0, 'other': 0, 'sd': 0, 'si': 0}
            
            # now look through all protein mutations in the fraction, expected and not
            for prot_mutation, prot_entry in all_references[background][fraction].items():
                if prot_entry['total'] < cutoff:
                    continue
                    # find all DNA entries with high enough counts
                    # same protein mutation can arise from multiple DNA mutations
                    
                for dna_mutation, dna_c in prot_entry['dna'].items():

                    if dna_c < cutoff:
                        continue

                    # sometimes a mutations will have a sequencing error / silent substitution elsewhere, which biases the count
                    dna_mut_type = classify_dna(dna_mutation)
                    prot_mut_type = classify_protein(prot_mutation)

                    if dna_mut_type == 'wt':
                        continue
                    elif dna_mut_type == 'b':
                        continue
                    
                    # it's possible that the DNA mutation looks valid but there is something off on the protein level
                    # skip those

                    if prot_mut_type is None:
                        print("None", prot_mut_type, prot_mutation, dna_mut_type, dna_mutation)
                        continue
                    elif (dna_mut_type in ['d3', 'd6', 'd9', 'i3', 'i6', 'i9']) and (prot_mut_type not in translate_variant_type(dna_mut_type)):
                        continue                       

                    try:
                        dna_count[background + '.' + fraction][dna_mut_type] += 1
                        dna_reads[background + '.' + fraction][dna_mut_type] += dna_c
                        distinct_mutations += 1
                        total_count += dna_c                     
                    except KeyError:
                        dna_count[background + '.' + fraction]['other'] += 1
                        dna_reads[background + '.' + fraction]['other'] += dna_c


    return dna_count, dna_reads

def get_full_protein_composition(all_references, cutoff=10):
    """
    This version does no filtering for correspondence between DNA and protein mutations, just checking what we see in the library.
    Probably want to use cutoff >1 for deletion libraries.
    """
    protein_count = {}
    protein_reads = {}
    
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            distinct_mutations = 0
            total_count = 0
            
            protein_count[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            '*': 0, 'f': 0, 'other': 0, 'b': 0}
            protein_reads[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            '*': 0, 'f': 0, 'other': 0, 'b': 0}
            
            for prot_mutation, entry in all_references[background][fraction].items():
                mut_total = entry['total']
                if entry['total'] < cutoff:
                    continue
                    
                # find all DNA entries with high enough counts
                # make sure the contributing DNA mutatios don't involve spurious extra substitutions
                prot_type = classify_protein(prot_mutation)

                try:
                    protein_count[background + '.' + fraction][prot_type] += 1
                    protein_reads[background + '.' + fraction][prot_type] += mut_total
                    distinct_mutations += 1
                    total_count += mut_total
                except KeyError:
                    # print(prot_type, mutation, mut_total)
                    protein_count[background + '.' + fraction]['other'] += 1
                    protein_reads[background + '.' + fraction]['other'] += mut_total

    return protein_count, protein_reads

def get_filtered_protein_composition(all_references, cutoff=5):

    """
    This version does filter for 'allowed' dna and protein mutations
    Probably want to use cutoff >1 for deletion libraries. Need cutoff >5 for subsitutions to avoid sequencing noise.
    """
    protein_count = {}
    protein_reads = {}
    allowed_dna = ['d3', 'd6', 'd9', 'i3', 'i6', 'i9', 's']
    allowed_protein = [ 'd', 'c-sd', 'c-dd', 'c-sdd', 'c-ddd','c-sddd',
                       'i1','c-si1', 'i2','c-si2', 'i3','c-si3',
                       '*','s', 'c-ss'] # stop codons c
    
    for background in all_references.keys():
        for fraction in all_references[background].keys():
            distinct_mutations = 0
            total_count = 0
            
            protein_count[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            '*': 0}
            protein_reads[background + '.' + fraction] = {'s': 0, 'd': 0, 'c-dd': 0, 'c-ddd': 0, 'i1': 0, 'i2': 0,
                            'i3': 0, 'c-ss': 0, 'c-sd': 0, 'c-sdd': 0, 'c-sddd': 0, 'c-si1': 0, 'c-si2': 0, 'c-si3': 0,
                            '*': 0}
            
            for prot_mutation, entry in all_references[background][fraction].items():
                mut_total = entry['total']
                if entry['total'] < cutoff:
                    continue
                # any valid protein muttion should be within the ones expected in the fraction
                # find all DNA entries with high enough counts
                # make sure the contributing DNA mutatios don't involve spurious extra substitutions
                prot_type = classify_protein(prot_mutation)
                if prot_type not in allowed_protein:
                    continue

                # make sure that the protein mutation arises from what we consider a valid DNA mutation
                # need special handling for STOPs
                valid_dna = False
                for dna_error, dna_c in entry['dna'].items():
                    dna_type = classify_dna(dna_error)
                    
                    if prot_type == '*':
                        if fraction in allowed_dna:
                            if dna_type == fraction:
                                valid_dna = True
                                break
                        else:
                            if dna_type in allowed_dna:
                                valid_dna = True
                    elif (dna_type in allowed_dna) and prot_type in translate_variant_type(dna_type):
                        valid_dna = True
                        break
                
                if valid_dna:
                    protein_count[background + '.' + fraction][prot_type] += 1
                    protein_reads[background + '.' + fraction][prot_type] += mut_total
                    distinct_mutations += 1
                    total_count += mut_total

    return protein_count, protein_reads