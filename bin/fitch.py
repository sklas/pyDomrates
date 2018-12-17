def fitch(node):
    """Infers the state of arrangements for all ancestral nodes in a given tree 
    by majority rule.  Also tries to determine single-step events that could 
    lead to the gain and loss of arrangements.
    This function expects a bifurcating tree as an ete2-object, in which the 
    leafs contain the feature 'arrangements' (a set of all arrangements present
    in the leaf).
    """
    leaves_to_root(node)
    root_to_leaves(node)
    ancestral_states(node)

def parental_state(left_child, right_child):
    """Calculate the parental state of an domain arrangement.  The states for
    this arrangement in the two children are given to this function as integers
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

def leaves_to_root(node):
    """For one node, collect all the arrangements that are present in the 
    children.  If a single arrangement is present in both children, set the 
    state for this arrangement to present (1).  If an arrangement is only 
    present in one of the children, set the state to unknown (0).  Because we 
    are only comparing arrangements that are present in at least one of the 
    children, getting an arrangement which is not present in both children is 
    impossible.  
    """
    if node.is_leaf():
        return dict((arr, 1) for arr in node.arrangements.keys())

    state = dict()
    cl, cr = node.get_children()
    state_l = leaves_to_root(cl)
    state_r = leaves_to_root(cr)
    arrangements_l = set(state_l.keys())
    arrangements_r = set(state_r.keys())
    for arrangement in arrangements_l.union(arrangements_r):
        state[arrangement] = parental_state(state_l.get(arrangement, -1), 
            state_r.get(arrangement, -1))

    node.add_features(arrangements=state)
    return node.arrangements

def root_to_leaves(node):
    """Determine the state for uncertain arrangements in a node according to 
    the state it has in the parents.  If a domain status is uncertain at the 
    root, set it to present (1).  In other nodes, determine the status of 
    uncertain arrangements by looking that particular arrangement up in the
    parent node.  If it is present in the parent, set it to present (1),
    otherwise delete it from that nodes repository.
    """
    if node.is_leaf():
        return
    elif node.is_root():
        # delete all arrangements which are definitely not there and set the
        # status for unsure ones to present
        node.arrangements = set(a for (a, c) in node.arrangements.iteritems() 
                if c >= 0)
    else:
        # delete all arrangements that are definitely not present at this node
        # and all arrangement with state 0, if the domain does not exist in
        # the parent.
        for a, c in node.arrangements.items(): # work on a copy
            if c == -1:
                del node.arrangements[a]
            elif c == 0 and a not in node.up.arrangements:
                del node.arrangements[a]
        # count numbers are now irrelevant
        node.arrangements = set(node.arrangements)

    cl, cr = node.get_children()
    root_to_leaves(cl)
    root_to_leaves(cr)

def ancestral_states(node):
    """Calculate which arrangements have been gained or lost on all nodes in 
    the tree.
    """
    d = set(node.arrangements)
    if node.is_root():
        p = set()
    else:
        p = set(node.up.arrangements)
    gained = d.difference(p)
    lost = p.difference(d)
    node.add_features(gained_arr=gained, lost_arr=lost)

    if node.is_leaf():
        return

    cl, cr = node.get_children()
    ancestral_states(cl)
    ancestral_states(cr)

if __name__ == '__main__':
    import ete2
    import dollo
    tree = ete2.Tree('(((A,B),C),D);')
    (tree&'A').add_features(domains=set('ABCDE'), arrangements={'A;B;C': 1, 'A': 1, 'D;E': 1})
    (tree&'B').add_features(domains=set('ABCDE'), arrangements={'A;B;C': 1, 'A;B': 1, 'D':1, 'E':1}) 
    (tree&'C').add_features(domains=set('ABCDEF'), arrangements={'A;B': 2, 'D;E': 1, 'C;F': 1, 'C':1})
    (tree&'D').add_features(domains=set('ABCDE'), arrangements={'A;B': 2, 'D;E': 1, 'C':1})
    dollo.dollo(tree)
    fitch(tree)
