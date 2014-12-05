#!/usr/bin/env python

import argparse, newick3, phylo3, math, os, random, shutil, subprocess, sys, time
import numpy as np
from dendropy import treesim

results_dir = os.path.expanduser("~/taxonjackknife/data_products")
temp_dir = os.path.expanduser("~/temp")

equal_rates_dir = results_dir + "/equal_rates_model"
random_rates_dir = results_dir + "/random_rates_model"
pectinate_random_dir = results_dir + "/pectinate_random_model"
balanced_random_dir = results_dir + "/balanced_random_model"
balanced_random_short_tips_dir = results_dir + "/balanced_random_short_tips_model"
TEST_DIR = results_dir + "/TEST"

score_file_column_labels = ["j_freq","j_ica","b_freq","b_ica","length","depth","in_true_tree"]

n_random_trees = 5 # currently 5
n_tips_per_tree = 10 #1000 # currently 1000
n_reps_taxon_jackknife = "10" #"200" # "200"
n_reps_bootstrap = "10" #"100" # "100"

# simple simulation parameters
birth_rate = 1.0
death_rate = birth_rate / 10

expected_tree_depth = math.log(n_tips_per_tree) / (birth_rate - death_rate)

n_threads = "30"
n_threads_bootstrap = "9"
max_bs_branches_to_process = 10000 / n_random_trees # test max 10,000 non-existent branches

_j = 1
def recur_balanced(p, max_d, max_t, _d=0):
    
    global _j
    
    if max_t < 2 or _d > max_d:
        return

    _d += 1
    n_max_t = [math.floor(max_t/2),]
    n_max_t.append(max_t - n_max_t[0])
    for i in range(2):
        c = phylo3.Node()
        p.add_child(c)
        if n_max_t[i] < 2:
            c.istip = True
            c.label = "T"+str(_j)
            _j += 1
        else:
            recur_balanced(c, max_d, n_max_t[i], _d)

def get_balanced_tree(branch_lengths_function):

    global _j
    _j = 1
    
    tree_depth = math.floor(math.log(n_tips_per_tree,2))

    tree = phylo3.Node()
    recur_balanced(tree, tree_depth, n_tips_per_tree)

    mean_br_length = expected_tree_depth / tree_depth
    return branch_lengths_function(tree, mean_br_length)

def get_pectinate_tree(branch_lengths_function):

    tree = phylo3.Node()
    n = tree
    for i in range(n_tips_per_tree - 1):
        c = phylo3.Node()
        t = phylo3.Node()
        t.istip = True
        t.label = str(i+1)
        n.add_child(c)
        n.add_child(t)
        n = c
        if i == n_tips_per_tree - 2:
            c.istip = True
            c.label = "T"+str(i+2)
    
    mean_br_length = expected_tree_depth / (n_tips_per_tree - 1)
    return branch_lengths_function(tree, mean_br_length)

def assign_branch_lengths_from_tips(tree, mean_br_length):

    # NOTE: assumes a bifurcating tree

    for n in tree.iternodes(order=phylo3.POSTORDER):

        if n.istip:
        
            # check if this is a cherry
            p = n.parent
            for s in p.children:
                if s != n:
                    if s.istip and hasattr(s, "height"):
                        # a cherry, other tip has an assigned length, copy it
                        n.length = s.length
                    else: # not a cherry or other tip has no length, assign a length
                        n.length = np.random.exponential(mean_br_length, 1)[0]
            n.height = 0 # tips are always height 0

        else: # not a tip
        
            # collect lengths and heights of children
            l = (n.children[0].length, n.children[1].length)
            h = (n.children[0].height, n.children[1].height)

            # normalize the lengths of this node's children so the
            # cumulative distance to the tips is equal on both sides
            if h[0] + l[0] > h[1] + l[1]:
                m = 0
                n.children[1].length += h[0] + l[0] - (h[1] + l[1])
            else:
                m = 1
                n.children[0].length += h[1] + l[1] - (h[0] + l[0])
                
            # assign length and height to this node
            n.length = np.random.exponential(mean_br_length/2, 1)[0]
            n.height = n.children[m].height + n.children[m].length

    return newick3.to_string(tree)    

def assign_branch_lengths_from_root(tree, mean_br_length):

    # get random draws, one more than we need for just the internal branches
    l = [x for x in np.random.exponential(mean_br_length, n_tips_per_tree)]

    # assign the internal branch lengths    
    tips_seen = 0
    for i, n in enumerate(tree.iternodes()):
        if not n.istip:
            n.length = l[i - tips_seen]
        else:
            tips_seen += 1

    # find the longest root-tip path (i.e. shallowest cherry)
    youngest_cherry = None
    tree_depth = 0
    for t in tree.leaves():
        p = t
        pathlen = 0
        while p.parent != None:
            pathlen += p.length
            p = p.parent

        if pathlen > tree_depth:
            tree_depth = pathlen
            youngest_cherry = t.parent

    # use the last branch length for the children of the shallowest node
    for m in youngest_cherry.children:
        m.length = (l[-1] / 10)
        
    tree_depth += youngest_cherry.children[0].length
    
    # assign tip lengths    
    for n in tree.iternodes():

        # get the current node depth
        p = n
        d = p.length
        while p.parent != None:
            p = p.parent
            d += p.length
        
        for c in n.children:
            if c.istip:
                c.length = tree_depth - d
        
    return newick3.to_string(tree)

def get_random_tree(branch_lengths_function):

    # note: the branch lengths function is unused in this tree generation method

    simt = treesim.birth_death(birth_rate=birth_rate, death_rate=death_rate, ntax=n_tips_per_tree)
    return simt.as_newick_string()

# for the equal rates model (simplest case)
def simulate_equal_rates(tree_label, tree_function, branch_lengths_function):

    transition_rate = birth_rate / 10
    gtr_equal = " ".join([str(transition_rate),]*6)
        
    indelible_control_file_text = """\
[TYPE] NUCLEOTIDE 2	//  nucleotide simulation using algorithm from method 2.
[MODEL]    gtr_equal
  [submodel]  GTR {model}
  [statefreq] 0.25 0.25 0.25 0.25 //  pi_T=0.25, pi_C=0.25, pi_A=0.25, pi_G=0.25
[TREE] tree {tree}
[PARTITIONS] part   [tree gtr_equal {aln_length}] // alignment length of {aln_length}
[EVOLVE]
  part 1 {tree_label} //  1 replicate generated from partition 'part' in file '{tree_label}.fas'
"""

    # repeat until we get an acceptable tree
    while True:

        # randomly generate a tree
        tree_string = tree_function(branch_lengths_function)

        # save the tree topology to a file
        with open(tree_label + ".tre", "w") as tree_file:
            tree_file.write(tree_string)
            tree_file.write(";")

        # write a control file for indelible
        with open("control.txt","w") as control_file:
            control_file.write(indelible_control_file_text.format(tree="("+tree_string+");", model=gtr_equal, aln_length=10000, tree_label=tree_label))

        # simulate data on the tree
        p = subprocess.Popen("indelible", stdout=subprocess.PIPE)
        r = p.communicate()
        if r[0].find("ERROR in [TREE] block tree in command [user]") < 0:
            # there was no error (substring index == -1) so move on 
            break
        
    # return the tree file and alignment file names 
    return (tree_label + ".tre", "%s_TRUE.phy" % tree_label)

# for the equal rates model (simplest case)
def simulate_random_rates(tree_label, tree_function, branch_lengths_function):
        
    indelible_control_file_text = """\
[TYPE] NUCLEOTIDE 2	//  nucleotide simulation using algorithm from method 2.

[MODEL]    gtr1
  [submodel]  GTR {models[0]}
  [statefreq] {statefreqs[0]}

[MODEL]    gtr2
  [submodel]  GTR {models[1]}
  [statefreq] {statefreqs[1]}

[MODEL]    gtr3
  [submodel]  GTR {models[2]}
  [statefreq] {statefreqs[2]}

[MODEL]    gtr4
  [submodel]  GTR {models[3]}
  [statefreq] {statefreqs[3]}

[MODEL]    gtr5
  [submodel]  GTR {models[4]}
  [statefreq] {statefreqs[4]}

[TREE] tree {tree}
[PARTITIONS] part1   [tree gtr1 {part_length}]
[PARTITIONS] part2   [tree gtr2 {part_length}]
[PARTITIONS] part3   [tree gtr3 {part_length}]
[PARTITIONS] part4   [tree gtr4 {part_length}]
[PARTITIONS] part5   [tree gtr5 {part_length}]

[EVOLVE]
  part1 1 {tree_label}_part_1 //  1 replicate generated from partition 'part1' in file '{tree_label}_part_1.fas'
  part2 1 {tree_label}_part_2 //  1 replicate generated from partition 'part2' in file '{tree_label}_part_2.fas'
  part3 1 {tree_label}_part_3 //  1 replicate generated from partition 'part3' in file '{tree_label}_part_3.fas'
  part4 1 {tree_label}_part_4 //  1 replicate generated from partition 'part4' in file '{tree_label}_part_4.fas'
  part5 1 {tree_label}_part_5 //  1 replicate generated from partition 'part5' in file '{tree_label}_part_5.fas'
"""

    # simulation parameters
    min_transition_rate = birth_rate / 100
    max_transition_rate = birth_rate / 10
    min_state_freq = 0.1
    max_state_freq = 0.3
    part_length = 2000
    n_parts = 5
    aln_length = part_length * n_parts

    # repeat until we get an acceptable tree
    while True:
        
        models = []
        statefreqs = []
        for j in range(n_parts):
            m = []
            t = []
            for k in range(6):
            
                # first generate a transition rate
                g = -1
                while g < min_transition_rate or g > max_transition_rate:
                    g = random.random()
                m.append(g)

            # perturb the model more
            scalar = random.randint(2,20)/float(10)
            m = [c*scalar for c in m]
            models.append(" ".join([str(v) for v in m]))
            
            for k in range(3):
                # now generate a state frequency
                g = -1
                while g < min_state_freq or g > max_state_freq:
                    g = random.random()
                t.append(g)

            # calculate the final value
            t.append(str(1 - sum(t)))
            statefreqs.append(" ".join([str(v) for v in t]))

        # randomly generate a tree
        tree_string = tree_function(branch_lengths_function)

        # save the tree topology to a file
        with open(tree_label + ".tre", "w") as tree_file:
            tree_file.write(tree_string)
            tree_file.write(";")

        # write a control file for indelible
        with open("control.txt","w") as control_file:
            control_file.write(indelible_control_file_text.format(tree="("+tree_string+");", models=models, \
                    statefreqs=statefreqs, part_length=aln_length/5, tree_label=tree_label))

        # simulate data on the tree
        p = subprocess.Popen("indelible", stdout=subprocess.PIPE)
        r = p.communicate()
        if r[0].find("ERROR in [TREE] block tree in command [user]") < 0:
            # there was no error (substring index == -1) so move on 
            break
    
    # combine the alignments and produce a partitions file
    aln = {}
    for j in range(n_parts):
        with open(tree_label+"_part_"+str(j+1)+"_TRUE.phy","r") as p:
            data = read_phylip(p)
        for name, seq in data.items():
            if name not in aln:
                aln[name] = ""
            aln[name] += seq

    alignment_file_name = "%s_combined_aln.phy" % tree_label
    with open(alignment_file_name,"w") as alignment_file:
        alignment_file.write(str(n_tips_per_tree) + " " + str(aln_length) + "\n")
        for name, seq in aln.items():
            alignment_file.write(name + " " + seq + "\n")
    
    partitions_file_name = "%s_combined_part.txt" % tree_label
    with open(partitions_file_name,"w") as partitions_file:
        for j in range(n_parts):
            partitions_file.write("DNA, p{j} = {begin}-{end}\n".format(j=j, begin=j*part_length+1, end=j*part_length+part_length))
            
    # return the tree file, alignment file, and partitions file names
    return (tree_label + ".tre", alignment_file_name, partitions_file_name)

def read_phylip(infile):
    aln = {}
    first_line = True
    for line in infile:
        
        if first_line == True:
            first_line = False
            continue
        
        parts = line.split(" ")
        if len(parts) < 2:
            continue
        
        if len(parts[0]) < 1:
            continue
        
        s = ""
        j = 1
        while parts[j] != "\n":
            s += parts[j]
            j += 1
        
        aln[parts[0].strip()] = s

    return aln

def subsample_beta(path_to_target_alignment, path_to_target_partfile):
    """beta"""
    cmd = 'subsample_alignment_beta.py '
    cmd += path_to_target_alignment 
    cmd += ' partfile=' + path_to_target_partfile if path_to_target_partfile is not None else ''
    subprocess.call(cmd, shell='True')
    return path_to_target_alignment.rsplit('.',1)[0]+'.subsampled.phy'

def run_single_tree(base_dir, simulation_function, tree_function, tree_number, branch_lengths_function=None, subsampling_function=None):

    base_dir += "_" + str(n_random_trees) + "_trees_of_" + str(n_tips_per_tree) + "_tips"
    
    if subsampling_function is not None:
        base_dir += "_SUBSAMPLED_" + subsampling_function.__doc__
    
    if not os.path.exists(base_dir):
        sys.exit("need to initialize the model dir first")

    final_scores_dir = base_dir + "/all_scores"
    all_scores_file_name = final_scores_dir + "/ALL.scores.csv"
    all_times_file_name = final_scores_dir + "/ALL.times.csv"
        
    os.chdir(base_dir)

    tree_label = "tree" + str(tree_number)
    if os.path.exists(tree_label):
        shutil.rmtree(tree_label) 
    working_tree_dir = os.getcwd() + "/" + tree_label
    os.mkdir(working_tree_dir)
    os.chdir(working_tree_dir)

    # run the tree and alignment simulation, get the results    
    simulation_results = simulation_function(tree_label, tree_function, branch_lengths_function)
    tree_file_name = simulation_results[0]
    alignment_file_name = working_tree_dir + "/" + simulation_results[1]
    partitions_file_name = working_tree_dir + "/" + simulation_results[2] if len(simulation_results) > 2 else None
    
    # subsample the alignment if necessary
    if subsampling_function is not None:
        alignment_file_name = subsampling_function(alignment_file_name, partitions_file_name)

    # calculate support for the original tree based on the data, support values will be stored in node_scores.csv
    subsample_args = ["subsample_edge_quartets.py",
        "-t", tree_file_name,
        "-n", alignment_file_name,
        "-#", n_reps_taxon_jackknife,
        "-T", n_threads,
        "-e", temp_dir,
        "-X", raxml_executable]

    if partitions_file_name is not None:
        subsample_args += ["-q", partitions_file_name]

    start = time.time()
    subprocess.call(subsample_args)
    jackknife_time = time.time() - start

    print("\ntime for taxon jackknife: %.2f seconds" % jackknife_time)
    with open(all_times_file_name,"a") as timefile:
        timefile.write("%s,%.2f," % (tree_number, jackknife_time))
    
    # get the node scores from the taxon jackknife
    node_scores = {}
    with open("node_scores.csv","r") as node_scores_file:
        first_line = True
        for line in node_scores_file:
            if first_line:
                first_line = False
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue
            node_scores[parts[0]] = {"j_freq":parts[1],"j_ica":parts[2],"in_true_tree":"TRUE"}

    # make a painted tree for reviewing the taxon jackknife results
    paint_args = ["paint_branches.py",
        "-t", "RESULT.labeled.tre",
        "-c", os.path.expanduser("~/scripts/figtree_color_palettes/blue_gray_red_1_to_neg_1.csv"),
        "-l", "ica",
        "-d", "node_scores.csv",
        "-n", tree_label]
    subprocess.call(paint_args)

    # calculate bootstrap values for the original tree based on the data
    raxml_bootstrap_args = [raxml_pthreads_executable,
        "-n", tree_label,
        "-s", alignment_file_name,
        "-#", n_reps_bootstrap,
        "-x", "123",
        "-p", "123",
        "-m", "GTRCAT",
        "-T", n_threads_bootstrap]
    
    if partitions_file_name is not None:
        subsample_args += ["-q", partitions_file_name]

    start = time.time()
    subprocess.call(raxml_bootstrap_args)
    bootstrap_time = time.time() - start

    print("\ntime for bootstrap: %.2f seconds" % bootstrap_time)
    with open(all_times_file_name,"a") as timefile:
        timefile.write("%.2f\n" % bootstrap_time)

    # call pxbp to generate bipart frequencies from bootstrap replicates
    pxbp_args = ["pxbp", "-t", "RAxML_bootstrap." + tree_label]
    p = subprocess.Popen(pxbp_args, stdout=subprocess.PIPE)
    results = p.communicate()[0]

    # parse pxbp results
    bootstrap_biparts = {}
    for line in results.split("\n"):
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 2:
            continue
        ingroup_labels = tuple(sorted(parts[0].split()))
        freq = None
        ica = None
        if parts[1].strip() == "1":
            freq = "1"
            ica = "1"
        else:
            freq = parts[2]
            ica = parts[4]
        bootstrap_biparts[ingroup_labels] = (freq, ica)

    # TODO: call mrbayes to generate a bayesian posterior distribution

    # load RESULT.labeled.tre
    labeled_tree = None
    true_tree_biparts = set()
    with open("RESULT.labeled.tre","r") as labeled_tree_file:
        labeled_tree = newick3.parse(labeled_tree_file.readline())

    # for each node in the labeled tre
    leaves = labeled_tree.leaves()
    for n in labeled_tree.iternodes():

        # skip the root and the tip nodes
        if n.parent == None or n.istip:
            continue
        
        # special case: if root position is adjacent to a tip
        if n.label not in node_scores and n.parent.parent == None:
            # in this case, the bipart is always defined--if you remove the root then the
            # bipart is the tip vs. the rest of the tree, so the taxon jackknife won't score
            # it (and its label won't show up in the node scores).
            continue

        true_tree_biparts.add(tuple(sorted([t.label for t in n.leaves()]))) 
        
        # record the branch length
        node_scores[n.label]["length"] = n.length
    
        # record the depth
        c = n
        d = 0
        while True:
            c = c.children[0]
            d += c.length
            if c.istip:
                break
        node_scores[n.label]["depth"] = d

        # see if this node's bipart is represented in the bootstraps, if so get the frequency and the ica score
        s = tuple(sorted([t.label for t in n.leaves()]))
        if s in bootstrap_biparts:
            node_scores[n.label]["b_freq"] = bootstrap_biparts[s][0]
            node_scores[n.label]["b_ica"] = bootstrap_biparts[s][1]
        else:
            node_scores[n.label]["b_freq"] = "0"
            node_scores[n.label]["b_ica"] = "NA"

    # process all bootstrap trees, recording info for all branches not in the true tree, so we can calc ica for them
    bootstrap_trees = []
    with open("RAxML_bootstrap." + tree_label, "r") as bootstrap_file:
        for line in bootstrap_file:
            bootstrap_trees.append(newick3.parse(line))

    # pull out all the taxon quartets defined by branches in all bootstrap trees
    bootstrap_quartets = set()
    leaves = labeled_tree.leaves()
    all_taxon_labels = set([l.label for l in leaves])
    for bootstrap_tree in bootstrap_trees:
        for node in bootstrap_tree.iternodes():

            # skip tips and the root node itself
            if node.istip or node.parent == None:
                continue

            # determine if this branch is represented in the true tree
            ingroup_labels = tuple(sorted([t.label for t in node.leaves()]))
            process_branch = False
            if ingroup_labels not in true_tree_biparts:
                remaining_labels = tuple(sorted(all_taxon_labels - set(ingroup_labels)))
                if remaining_labels not in true_tree_biparts:
                    process_branch = True
        
            # if this branch isn't in the true tree, make a tree we can use to evaluate ica for it as well
            if process_branch:
                ########### modified from subsample_edge_quartets
        
                # get leaf sets for the four connected subtrees
            
                # two daughter subtrees
                r1 = set([node.children[0].label,] if node.istip else [l.label for l in node.children[0].leaves()])
                r2 = set([node.children[1].label,] if node.istip else [l.label for l in node.children[1].leaves()])

                # sibling/parent subtrees
                is_other_side_of_root = False # used when we hit the root for the second time
                skip_tip_child_of_root = False # used when one of the children of the root node is a tip
                tip_child_label = None
                for sib in node.parent.children:
                    if sib != node:

                        # if one of the subtrees is the root, skip over it
                        if len(sib.leaves()) + len(node.leaves()) == len(leaves):

                            # if we already processed this bipart (on other side of the root), don't do it again
                            if (root_bipart_label != None):
                                is_other_side_of_root = True
                                break

                            # get the subtrees opposite the root
                            if len(sib.children) == 2:
                                l1 = set([sib.children[0].label,] if sib.children[0].istip else [l.label for l in sib.children[0].leaves()])
                                l2 = set([sib.children[1].label,] if sib.children[1].istip else [l.label for l in sib.children[1].leaves()])
                            elif len(sib.children) == 0:
                                skip_tip_child_of_root = True
                                tip_child_label = sib.label
                            else:
                                print("Node %s does not have exactly 2 children. It will be skipped." % k)
                                continue

                            # remember that we've already done the root, so we can skip it when we hit the other side
                            root_bipart_label = node.label

                        # otherwise not at root, all connected subtrees have children
                        else:

                            # sibling subtree
                            l1 = set([l.label for l in sib.leaves()])

                            # the rest of the tree
                            l2 = set()
                            for label in [l.label for l in leaves]:
                                if label not in r1 and \
                                   label not in r2 and \
                                   label not in l1:
                                        l2.add(label)
                    
                if skip_tip_child_of_root:
                    print("not evaluating tip child '" + tip_child_label + "' of the root (ica is 1.0, as for all tips).")
                    continue

                ######## end modified from subsample_edge_quartets.py
            
                q = (frozenset((tuple(sorted(l1)),tuple(sorted(l2)))),frozenset((tuple(sorted(r1)),tuple(sorted(r2)))))
            
                # make sure this quartet isn't already in the set in reverse before adding it
                if (q[1],q[0]) not in bootstrap_quartets:
                    bootstrap_quartets.add(q)

    print("Found %s bootstrap quartets not in the original tree" % len(bootstrap_quartets))

    # get ready to process the bad branch quartets
    missing_topos_dir = temp_dir + "/" + tree_label + "_bootstrap_quartets_missing_from_true_tree"
    if os.path.exists(missing_topos_dir):
        shutil.rmtree(missing_topos_dir)
    os.mkdir(missing_topos_dir)
    os.chdir(missing_topos_dir)

    # set the sampling interval so we only process max_bs_branches_to_process, evenly distributed across the input set
    sample_freq = len(bootstrap_quartets)/float(max_bs_branches_to_process)
    max_br = min(len(bootstrap_quartets),max_bs_branches_to_process)
    cur_br = 1
        
    for p, q in enumerate(bootstrap_quartets): # for p in range(x)

        # skip branches that don't fall within our sample based on max_bs_branches_to_process
        if not p % sample_freq < 1:
            continue

        print("on bootstrap quartet %s out of %s" % (cur_br, max_br))
        bad_branch_label = str(p) + "b"
        os.mkdir(bad_branch_label)
        os.chdir(bad_branch_label)

        cur_br += 1
            
        tip_label_sets = []
        with open(bad_branch_label+".tre","w") as topo_file:
            clades = []
            for j in q:
                for k in j:
                    tip_label_sets.append(k) 
            for ls in tip_label_sets:
                text = ",".join(ls)
                if len(ls) > 1:
                    text = "(" + text + ")"
                clades.append(text)
            newick = "((" + clades[0] + "," + clades[1] + "),(" + clades[2] + "," + clades[3] + "));"
            topo_file.write(newick)

        # calculate support for the bootstrap-inferred FALSE branch, support values will be stored in node_scores.csv
        subsample_args = ["subsample_edge_quartets.py",
            "-t", bad_branch_label+".tre",
            "-n", alignment_file_name,
            "-#", n_reps_taxon_jackknife,
            "-T", n_threads,
            "-e", temp_dir,
            "-p", "1", # specify a stop node number that we only want to process the first (i.e. the root) node
            "-X", raxml_executable]
        subprocess.call(subsample_args)     

        # extract the ica score for this node. we only care about the first bipart, so we read the second line in the file
        with open("node_scores.csv","r") as node_scores_file:
            parts = [t.strip() for t in node_scores_file.readlines()[1].split(",")]
            node_scores[bad_branch_label] = {"j_freq":parts[1],"j_ica":parts[2]}

        try:
            ingroup_labels = tuple(sorted(tip_label_sets[0] + tip_label_sets[1]))
            node_scores[bad_branch_label]["b_freq"] = bootstrap_biparts[ingroup_labels][0]
        except KeyError:
            ingroup_labels = tuple(sorted(tip_label_sets[2] + tip_label_sets[3]))
            node_scores[bad_branch_label]["b_freq"] = bootstrap_biparts[ingroup_labels][0]

        node_scores[bad_branch_label]["b_ica"] = bootstrap_biparts[ingroup_labels][1]
        node_scores[bad_branch_label]["in_true_tree"] = "FALSE"
        node_scores[bad_branch_label]["length"] = "NA"
        node_scores[bad_branch_label]["depth"] = "NA"

        os.chdir(missing_topos_dir)
        
    os.chdir(working_tree_dir)

    # write scores to file and prepare for next iteration
    with open(all_scores_file_name,"a") as scores_file:
        for node_label, values in node_scores.items():
            scores_file.write(str(tree_number)+","+node_label)
            for c in score_file_column_labels:
                scores_file.write(","+str(values[c]))
            scores_file.write("\n")

def init_model(base_dir, subsampling_function=None):
    
    base_dir += "_" + str(n_random_trees) + "_trees_of_" + str(n_tips_per_tree) + "_tips"
    if subsampling_function is not None:
        base_dir = base_dir + "_SUBSAMPLED_" + subsampling_function.__doc__

    x = raw_input("Are you sure you want to initialize '" + base_dir + "'? y/n: ")
    if x[0].lower() != "y":
        exit(0)

    final_scores_dir = base_dir + "/all_scores"
    all_scores_file_name = final_scores_dir + "/ALL.scores.csv"
    all_times_file_name = final_scores_dir + "/ALL.times.csv"

    for d in [base_dir, final_scores_dir]:
        if not os.path.exists(d):
            os.mkdir(d)
    os.chdir(base_dir)

    with open(all_scores_file_name,"w") as scores_file:
        scores_file.write("tree,node,"+",".join(score_file_column_labels)+"\n")

    with open(all_times_file_name,"w") as times_file:
        times_file.write("tree,jackknife_time,bootstrap_time\n")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="taxon jackknife simulations")
    
    parser.add_argument("-s", "--raxml-single", nargs=1, required=True, help="The name of the raxmlHPC single thread executable.")

    parser.add_argument("-p", "--raxml-pthreads", nargs=1, required=True, help="The name of the raxmlHPC PTHREADS executable.")

    args = parser.parse_args()

    raxml_executable = args.raxml_single[0]
    raxml_pthreads_executable = args.ramlx_pthreads[0]

    if len(sys.argv) > 1:
        raxml_executable = sys.argv[1]
        raxml_pthreads_executable = sys.argv[2]
    
    for d in [temp_dir, results_dir]:
        if not os.path.exists(d):
            os.mkdir(d)

    sys.setrecursionlimit(1500)

    # A    
#    init_model(equal_rates_dir)
#    for i in range(1, n_random_trees+1):
#       run_single_tree(equal_rates_dir, simulate_equal_rates, get_random_tree, i)

#    # B
#    init_model(random_rates_dir)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(random_rates_dir, simulate_random_rates, get_random_tree, i)

    # C
#    init_model(pectinate_random_dir)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(pectinate_random_dir, simulate_random_rates, get_pectinate_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_root)

    # D
#    init_model(balanced_random_dir)
#    for i in range(3, n_random_trees+1):
#        run_single_tree(balanced_random_dir, simulate_random_rates, get_balanced_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_root)

    # E
#    init_model(balanced_random_short_tips_dir)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(balanced_random_short_tips_dir, simulate_random_rates, get_balanced_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_tips)

    # test
#    init_model(TEST_DIR, subsample_beta)
#    for i in range(1, 20):
#        run_single_tree(TEST_DIR, simulate_random_rates, get_balanced_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_tips, \
#                        subsampling_function=subsample_beta)

### using beta subsampling

    # A    
#    init_model(equal_rates_dir, subsampling_function=subsample_beta)
#    for i in range(1, n_random_trees+1):
#       run_single_tree(equal_rates_dir, simulate_equal_rates, get_random_tree, i, subsampling_function=subsample_beta)

#    # B
#    init_model(random_rates_dir, subsampling_function=subsample_beta)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(random_rates_dir, simulate_random_rates, get_random_tree, i, subsampling_function=subsample_beta)

    # C
#    init_model(pectinate_random_dir, subsampling_function=subsample_beta)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(pectinate_random_dir, simulate_random_rates, get_pectinate_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_root, subsampling_function=subsample_beta)

    # D
#    init_model(balanced_random_dir, subsampling_function=subsample_beta)
#    for i in range(3, n_random_trees+1):
#        run_single_tree(balanced_random_dir, simulate_random_rates, get_balanced_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_root, subsampling_function=subsample_beta)

    # E
#    init_model(balanced_random_short_tips_dir, subsampling_function=subsample_beta)
#    for i in range(1, n_random_trees+1):
#        run_single_tree(balanced_random_short_tips_dir, simulate_random_rates, get_balanced_tree, i, \
#                        branch_lengths_function=assign_branch_lengths_from_tips, subsampling_function=subsample_beta)

    # test
    init_model(TEST_DIR, subsampling_function=subsample_beta)
    for i in range(1, 20):
        run_single_tree(TEST_DIR, simulate_random_rates, get_balanced_tree, i, \
                        branch_lengths_function=assign_branch_lengths_from_tips, \
                        subsampling_function=subsample_beta)
