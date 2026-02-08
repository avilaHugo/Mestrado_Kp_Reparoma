import Bio
from Bio import Phylo
from Bio.Phylo.NewickIO import _translate_name
read in newick file

tree = Phylo.read("tree.nwk", "newick")
open txt file with name replacements

with open("name_replacements.txt") as f:
replacements = {}
for line in f:
original, replacement = line.strip().split("\t")
replacements[original] = replacement
translate names in newick tree using name replacements

for clade in tree.find_clades():
if clade.name:
clade.name = _translate_name(replacements, clade.name)
write modified tree to new file

Phylo.write(tree, "modified_tree.nwk", "newick")

print("Name replacements complete!")

