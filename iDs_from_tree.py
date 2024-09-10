from ete3 import Tree

# Replace 'final_merged.ph' with the path to your .ph file
tree_file = 'final_merged.ph'

# Read the tree
tree = Tree(tree_file, format=1)

# Extract seqIDs that start with "k101"
k101_seqIDs = [leaf.name for leaf in tree.iter_leaves() if leaf.name.startswith("k101")]

# save the k101 seqIDs to a file
with open('extracted_seqIDs_from_ph_file.txt', 'w') as f:
    for seqID in k101_seqIDs:
        f.write(f"{seqID}\n")


