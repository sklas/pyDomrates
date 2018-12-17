# pyDomrates
pyDomrates is the python implemented prototype of domrates, a tool to reconstruct domain and domain arrangement content at the inner nodes of a phylogenetic tree.
With pyDomrates you can infer four main types of events being responsible for new domain arrangements, i.e.:
1. fusion of two arrangemnts to one novel,
2. fission of one arrangement into two novel,
3. loss of domains in an arrangement and
4. addition of novel domains.

Details of the algorithm can be found in (cite moore, kersting, domrates(submitted, maybe biorxive?)).

**We strongly encourage you to use our new C++ implementation of domrates! Check out https://domainworld.uni-muenster.de/programs/domrates/ and https://ebbgit.uni-muenster.de/domainWorld/DomRates**

The tool takes three inputs: 
* A phylogenetic tree of the species of interest,
* the path to the domain annotation for each species and
* the name of two leafs of the tree to determine the outgroup.

Domain annotations are the output of the pfamScan.pl script provided by Pfam.
pyDomrates outputs the domain and arrangement content for each node, the number of events used and  an annotated pylogenetic tree.

A typical call to the program might look like:
``` bash
PYTHONPATH=/home/user/domrates/lib python /home/user/domrates/domaintree.py -d /home/user/project/pfam -t /home/user/project/tree.nwk -c "Daphnia pulex" "Apis mellifera" -e ".pfam"
```
