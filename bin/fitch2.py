def fitch2(node):
    """Infers the state of domains for all ancestral nodes in a given tree 
    by majority rule.  Also tries to determine single-step events that could 
    lead to the gain and loss of domains.
    This function expects a bifurcating tree as an ete2-object, in which the 
    leafs contain the feature 'domains' (a set of all domains present
    in the leaf).
    """
    leaves_to_root2(node)
    root_to_leaves2(node)
    ancestral_states2(node)

def parental_state2(left_child, right_child):
    """Calculate the parental state of an domain domain.  The states for
    this domain in the two children are given to this function as integers
    taking only the values 1 (+, certainly there), 0 (?, status unclear) and -1
    (-, definitely not there).  Arithmetics look like this:
    +, + => +        1 +  1 =  1  <<
    +, ? => +        1 +  0 =  1
    ?, + => +        0 +  1 =  1
    -, - => -       -1 + -1 = -1  <<
    -, ? => -       -1 +  0 = -1
    ?, - => -        0 + -1 = -1
    -, + => ?       -1 +  1 =  0
    +, - => ?        1 + -1 =  0
    ?, ? => ?        0 +  0 =  0  <<

    This function is used to infer the ancestral states on the way from leaves
    to root.
    """
    if left_child == right_child:
        return left_child
    else:
        return left_child + right_child

def leaves_to_root2(node):
    """For one node, collect all the domains that are present in the 
    children.  If a single domain is present in both children, set the 
    state for this domain to present (1).  If an domain is only 
    present in one of the children, set the state to unknown (0).  Because we 
    are only comparing domains that are present in at least one of the 
    children, getting an domain which is not present in both children is 
    impossible.  
    """
    if node.is_leaf():
        return dict((arr, 1) for arr in node.domains)

    state = dict()
    cl, cr = node.get_children()
    state_l = leaves_to_root2(cl)
    state_r = leaves_to_root2(cr)
    domains_l = set(state_l)
    domains_r = set(state_r)
    for domain in domains_l.union(domains_r):
        state[domain] = parental_state2(state_l.get(domain, -1), 
            state_r.get(domain, -1))

    node.add_features(domains=state)
    return node.domains

def root_to_leaves2(node):
    """Determine the state for uncertain domains in a node according to 
    the state it has in the parents.  If a domain status is uncertain at the 
    root, set it to present (1).  In other nodes, determine the status of 
    uncertain domains by looking that particular domain up in the
    parent node.  If it is present in the parent, set it to present (1),
    otherwise delete it from that nodes repository.
    """
    if node.is_leaf():
        return
    elif node.is_root():
        # delete all domains which are definitely not there and set the
        # status for unsure ones to present
        node.domains = set(a for (a, c) in node.domains.iteritems() 
                if c >= 0)
    else:
        # delete all domains that are definitely not present at this node
        # and all domain with state 0, if the domain does not exist in
        # the parent.
        for a, c in node.domains.items(): # work on a copy
            if c == -1:
                del node.domains[a]
            elif c == 0 and a not in node.up.domains:
                del node.domains[a]
        # count numbers are now irrelevant
        node.domains = set(node.domains)

    cl, cr = node.get_children()
    root_to_leaves2(cl)
    root_to_leaves2(cr)

def ancestral_states2(node):
    """Calculate which domains have been gained or lost on all nodes in 
    the tree.
    """
    d = set(node.domains)
    if node.is_root():
        p = set()
    else:
        p = set(node.up.domains)
    gained = d.difference(p)
    lost = p.difference(d)
    node.add_features(gained_domains=gained, lost_domains=lost)

    if node.is_leaf():
        return

    cl, cr = node.get_children()
    ancestral_states2(cl)
    ancestral_states2(cr)

if __name__ == '__main__':
    import ete2
    import dollo
    tree = ete2.Tree('(((A,B),C),D);')
    #(tree&'A').add_features(domains=set('ABCDE'), domains={'A;B;C': 1, 'A': 1, 'D;E': 1})
    #(tree&'B').add_features(domains=set('ABCDE'), domains={'A;B;C': 1, 'A;B': 1, 'D':1, 'E':1}) 
    #(tree&'C').add_features(domains=set('ABCDEF'), domains={'A;B': 2, 'D;E': 1, 'C;F': 1, 'C':1})
    #(tree&'D').add_features(domains=set('ABCDE'), domains={'A;B': 2, 'D;E': 1, 'C':1})
    dollo.dollo(tree)
    fitch(tree)

