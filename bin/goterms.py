import os
# import collections
# import rpy2.robjects as robjects
# import sys
# import pytagcloud
# import math
'''
make an enrichment analysis of gained/lost genes with
topGO and rpy2
'''

"""make domain2go mappings and assign the universe. The universe 
includes all GO terms present in all domains in the tree, 
specified by genes variable. Output are two dicts, first with 
domain to GO term mapping, second with GO term to domain mapping
"""
def make_domain_go_mappings(mapping_file, genes, db):
    gene_to_go = collections.defaultdict(list)
    go_to_gene = collections.defaultdict(list)
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.startswith("!"):
                continue
            line = line.split()
            if db == "pfam":
                gene_id = line[0].split(":")[1]
                go_id = line[-1].strip()
            elif db == "supfam":
                gene_id = str(line[1])
                go_id = line[2]
            if gene_id in genes:
                gene_to_go[gene_id].append(go_id)
                go_to_gene[go_id].append(gene_id)
    return dict(gene_to_go), dict(go_to_gene)


"""make topGO object and run an analysis with rpy2.
Arguments are the domain to GO term mapping (gene_to_go), 
a list of domains of interest (subset), e.g. lost or gained 
domains, the type of go terms (e.g. "MF" for Molecular function), 
the p-value for significance and a boolean value for printing the results.
"""
def make_topGO_object(gene_to_go, subset, go_term_type, go_pval, printtable):
    # temporarily redirect standard output to nowhere, keeps the output file clean
    class Blackhole(object):
        def write(self, string):
            pass
    stdout = sys.stdout
    sys.stdout = Blackhole()
    def _set_to_namedvector(init_set):
        return robjects.r.c(*init_set)
    def _dict_to_namedvector(init_dict):
        return robjects.r.c(**init_dict)
    with open('tmp/rgenelist','w') as h:
        for gene, gos in gene_to_go.iteritems():
            h.write(gene+"\t"+",".join(gos)+"\n")

    # domains_in_universe = robjects.r.names(rgene_to_go)
    # rdoi = _set_to_namedvector(doi)
    robjects.r("library(topGO)")
    robjects.globalenv["subset"] = _set_to_namedvector(subset)
    #robjects.globalenv["domains_in_universe"] = domains_in_universe
    #robjects.globalenv["gene_to_go"] = _dict_to_namedvector(gene_to_go)
    robjects.globalenv["gotermtype"] = go_term_type
    #for dom in subset:
    #    if dom in gene_to_go:
    #       print dom
    #    else:
    #        print "???????"
    #print "\n\n\n"
    go_data = robjects.r('''
        GOmap <- readMappings(file="tmp/rgenelist")
        domains_in_universe <- names(GOmap)
        doi <- factor(as.integer(domains_in_universe %in% subset))
        names(doi) <- domains_in_universe
        new("topGOdata", ontology = gotermtype, allGenes = doi, annot = annFUN.gene2GO, gene2GO = GOmap)
    ''')
    os.remove("tmp/rgenelist")
    results = robjects.r.runTest(go_data, algorithm = "weight01", statistic = "fisher")
    scores = robjects.r.score(results)
    score_names = robjects.r.names(scores)
    num_summarize = min(50, len(score_names))
    results_table = robjects.r.GenTable(go_data, weight=results, orderBy="weight", topNodes=num_summarize)
   
    # standard out to stdout
    sys.stdout = stdout

    if printtable:
        print results_table
    # extract term names from the topGO summary dataframe
    GO_ID_INDEX = 0
    TERM_INDEX = 1
    ids_to_terms = dict()
    for index, go_id in enumerate(results_table[GO_ID_INDEX]):
        ids_to_terms[go_id] = results_table[TERM_INDEX][index]
    go_terms = []
    # convert the scores and results information into terms to return
    for index, item in enumerate(scores):
        if item < go_pval:
            go_id = score_names[index]
            go_terms.append((item, go_id, ids_to_terms.get(go_id, "")))
        go_terms.sort()
    return go_terms

def print_tag_cloud(go_term_list, nodename, db):
    tag_list = []
    for term in go_term_list:
        tag_list.append((term[2], (-math.log10(term[0]))))
    tags = pytagcloud.make_tags(tag_list, minsize=5, maxsize=30)
    pytagcloud.create_tag_image(tags, "results/tagclouds/"+db+"tagcloud_"+nodename+".png", size=(1000,800), layout=0)

def enrich_with_topGO(domains_of_interest, gene2go, go2gene, printtable=False):
    go_pval =  0.05 
    go_term_type = "BP"

    doi = set()
    
    go_terms = make_topGO_object(gene2go, domains_of_interest, go_term_type, go_pval, printtable)
    return go_terms
