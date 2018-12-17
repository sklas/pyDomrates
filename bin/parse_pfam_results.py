#!/usr/bin/env python2

import os

def load_res_file(filename, pfam_mapping=False,check_for_duplicates=False):
    """loads a pfamscan .res output file and extracts all domains that have
    an evalue less than a given threshold.  It extracts some information and
    returns a dict with protein ids as keys, the values are lists of tuples
    containing the domain name, alignment start/end, clan, evalue and domain
    type."""

    EVALUE_THRESHOLD = 1
    protein_hits = dict()

    with open(filename, 'r') as f:
        for line in f:
            # strip header
            if line[0] == '#' or not line.strip():
                continue

            elements = line.split()
            protein_id = elements[0]
            align_start = int(elements[1])
            align_end = int(elements[2])
            hmm_acc = elements[5].split('.')[0] # strip hmm family name, e.g. PF00105.13 -> PF00105
            hmm_name = elements[6]
            domain_type = elements[7]
            e_value = float(elements[12])
            clan = elements[14]
            if pfam_mapping:
                pfam_mapping[hmm_acc] = hmm_name

            if float(e_value > EVALUE_THRESHOLD):
                continue
            data = (hmm_acc, align_start, align_end, domain_type, clan)
            if protein_id not in protein_hits:
                protein_hits[protein_id] = [data]
            else:
                protein_hits[protein_id].append(data)
        
        # sort domains inside a protein according to sequence
        for protein, domains in protein_hits.items(): # work on a copy
            protein_hits[protein] = sorted(protein_hits[protein], key=lambda rank: (rank[1], rank[2]))

        # check for overlaps in annotation
        # TODO: too much overlaps, filter them out!
        if check_for_duplicates:
            f = os.path.basename(filename).split('-')[0]
            for protein, domains in protein_hits.iteritems():
                last_start = -1
                last_stop = -1
                for domain in domains:
                    start = domain[1]
                    stop = domain[2]
                    if start < last_stop:
                        # overlap
                        if stop < last_stop:
                            print('{}:{}: inclusion of domain {}'.format(f, protein, domain[0]))
                        else:
                            print('{}:{}: overlap of domain {}, {} residues'.format(f, protein, domain[0], last_stop-start+1))
                    last_start = start
                    last_stop = stop

        return protein_hits


def get_domain_data(proteins, clan=False, collapse=True):
    """Collect and return all domains and all domain arrangements present in 
    one organism.  For architectures, also count the occurences.
    Options include replacing domain names with pfam clan names, if available,
    and collapsing subsequent equal domains into one (only for arrangements).
    """

    domain_set = set()
    domain_counts = dict()
    genes_with_domain = 0
    architectures = dict()
    for protein, domains in proteins.iteritems():
        architecture = []
        last_domain = ''
        for domain in domains:
            domain_name = domain[0]
            clan_name = domain[4]
            if clan and clan_name.lower() != 'No_clan'.lower():
                domain_name = clan_name

            domain_set.add(domain_name)
            ## Uncomment for repeat domain counts on protein
            #domain_counts[domain_name] = domain_counts.get(domain_name, 0) + 1
            # create domain arrangement
            if collapse and domain_name == last_domain:
                continue
            domain_counts[domain_name] = domain_counts.get(domain_name, 0) + 1

            architecture.append(domain_name)
            last_domain = domain_name
        # collect and count architectures for organism
        architecture = ';'.join(architecture)
        architectures[architecture] = architectures.get(architecture, 0) + 1
    return (domain_set, architectures, domain_counts)
