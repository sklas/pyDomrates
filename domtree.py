#!/usr/bin/env python2

# Calculate domain content on all nodes in i5k tree

import ete3
import os
import sys
import parse_pfam_results as pfam
# import dollo
import fitch
import fitch2
from pylab import *
import argparse
from scipy import stats as stats
from collections import defaultdict
import numpy as np

def run(annotation_data, tree_filename, children, ending):
    # declare variables
    proteomes = {}

    # read trees
    f = open(tree_filename, 'r')
    tree_full_names = ete3.Tree(f.read(), format=1)
    f.close()

    # read data from pfam
    data_path = annotation_data
    species_list = dict()
    pfam_mapping = dict()
    for f in os.listdir(data_path):
        if f.endswith(ending):
            species_name = f.split('.')[0]
            content = pfam.load_res_file(os.path.join(data_path, f), pfam_mapping)
            proteomes[species_name] = content

    # build tree for the data
    tree = tree_full_names.copy()
    # assign data to the leafs
    all_domains_in_tree = set()
    for species_name, species_proteome in proteomes.iteritems():
        domains, arrangements, domain_counts = pfam.get_domain_data(species_proteome)
        all_domains_in_tree.update(domains)
        (tree&species_name).add_features(domains=domains, arrangements=arrangements, domain_counts=domain_counts)

    # Apply fitch parsimony
    fitch.fitch(tree)
    fitch2.fitch2(tree)

    # get domain wise counts
    dom_counts = {}
    for node in tree.iter_leaves():
        dom_counts[node.name] = node.domain_counts
    # make domainwise contingency table
    cont_table = defaultdict(dict)
    for species, counts in dom_counts.iteritems():
        for domain, count in counts.iteritems():
            cont_table[domain][species] = count

    # for species_name, proteome in proteomes.iteritems():
    #     print species_name
    #     print len(proteome)
    

    
    print "Direction, Age, Domain, SpeciesOfInterest, SpeciesComp, SOIIn, SOIOut, SCIn, SCOut, P-Value"
    for soi in ["Athalia_rosae", "Orussus_abietinus"]:
        for domain, acounts in cont_table.iteritems():
            for age in ["old", "young"]:
                if soi in acounts:
                    counts = acounts.copy()
                    ocounts = {}
                    if "Tribolium_castaneum" in counts:
                        ocounts["Tribolium_castaneum"] = counts.pop("Tribolium_castaneum")
                    if "Zootermopsis_nevadensis" in counts:
                        ocounts["Zootermopsis_nevadensis"] =counts.pop("Zootermopsis_nevadensis")
                    if soi == "Orussus_abietinus":
                        if "Athalia_rosae" in counts:
                            #ocounts["Athalia_rosae"] = 
                            counts.pop("Athalia_rosae")
                    elif soi == "Athalia_rosae":
                        if "Orussus_abietinus" in counts:
                            counts.pop("Orussus_abietinus")
                    if age == "old":
                        ocounts[soi] = counts[soi]
                        counts = ocounts
                    counts_sort = sorted(counts.values())
                    if not len(counts_sort) == 1:
                        if counts[soi] == counts_sort[0] and counts[soi] != counts_sort[1]:
                            for species, count in counts.iteritems():
                                if counts_sort[1] == count:
                                    chi2_data = np.array([[counts[soi], len(proteomes[soi])-counts[soi]], [count, len(proteomes[species])-count]])
                                    res = stats.chi2_contingency(chi2_data)
                                    # print res[1]
                                    if res[1] < 0.05:
                                        print ", ".join(["low", age, domain, soi, species, str(chi2_data[0,0]), str(chi2_data[0,1]), str(chi2_data[1,0]), str(chi2_data[1,1]), str(res[1])])
                                        break
                        if counts[soi] == counts_sort[-1] and counts[soi] != counts_sort[-2]:
                            for species, count in counts.iteritems():
                                if counts_sort[-2] == count:
                                    chi2_data = np.array([[counts[soi], len(proteomes[soi])-counts[soi]], [count, len(proteomes[species])-count]])
                                    res = stats.chi2_contingency(chi2_data)
                                    if res[1] < 0.05:
                                        print ", ".join(["high", age, domain, soi, species, str(chi2_data[0,0]), str(chi2_data[0,1]), str(chi2_data[1,0]), str(chi2_data[1,1]), str(res[1])])
                                        break

## Print tree
def print_tree(tree):
    # set node style
    nstyle = ete3.NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 0
    nstyle["fgcolor"] = "darkred"
    nstyle["vt_line_width"] = 7
    nstyle["hz_line_width"] = 7

    for node in tree.traverse():
        if node.is_leaf():
            # node.add_face(ete3.TextFace(len(node.arrangements), fsize=30), column=1, position="aligned")
            node.add_face(ete3.TextFace(len(node.domains), fsize=30), column=1, position="aligned")
        else:
            # node.add_face(ete3.TextFace(str(len(node.arrangements))+" ", fsize=25), column=0, position="branch-bottom")
            node.add_face(ete3.TextFace(str(len(node.domains))+" ", fsize=25), column=0, position="branch-bottom")
            node.add_face(ete3.TextFace(node.name+" ", fsize=25), column=0, position="branch-bottom")
        node.add_face(ete3.TextFace("+"+str(len(node.gained_domains)), fgcolor="green", fsize=25), column=0, position="branch-top")
        # node.add_face(ete3.TextFace("+"+str(len(node.gained_arr)), fgcolor="green", fsize=25), column=0, position="branch-top")
        node.add_face(ete3.TextFace("-"+str(len(node.lost_domains)), fgcolor="red", fsize=25), column=0, position="branch-top")
        # node.add_face(ete3.TextFace("-"+str(len(node.lost_arr)), fgcolor="red", fsize=25), column=0, position="branch-top")
        node.set_style(nstyle)

    def layout(node):
        N = ete3.AttrFace("name", fsize=30)
        N.margin_right = 10
        if node.is_leaf():
            if node.name == "Orussus_abietinus":
                N.background.color = "lightblue"
            elif node.name == "Athalia_rosae":
                N.background.color = "lightgreen"

            ete3.faces.add_face_to_node(N, node, 0, position="aligned")

    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.draw_guiding_lines = True
    # ts.guiding_lines_type = 0
    ts.extra_branch_line_color = "black"

    ts.extra_branch_line_type = 1
    ts.show_scale = False
    ts.layout_fn = layout
    ts.optimal_scale_level = "mid"

    tree.render("tree_gain_loss_domains_20160211.pdf", tree_style=ts, w=199, units="mm")


# Print Set 
def print_set(tree):
    all_doms = set()
    for node in tree.traverse():
        if node.name not in ("NoName" or ""):
            print(node.name)
            with open("domains/loss/"+str(node.name)+".doms", 'w') as f:
                f.write(("\n").join(node.lost_domains))
            all_doms |= node.domains


        print(node.name)
        with open("domains/all", 'w') as f:
            f.write(("\n").join(all_doms))


def main():
    """
        Compares trees with each other
    """
    parser = argparse.ArgumentParser(description='This is a program to assign domain content to species trees')
    parser.add_argument('-d', '--domainDir', help="Directory with domain files")
    parser.add_argument('-t', '--tree', help="The tree to use")
    parser.add_argument('-c', '--children', nargs=2, help='Children to use for root determination')
    parser.add_argument('-e', '--ending', default='.dom', help='The domain file ending')
    args = parser.parse_args()
    
    run(args.domainDir, args.tree, args.children, args.ending) 
    

if __name__ == "__main__":
    main()
