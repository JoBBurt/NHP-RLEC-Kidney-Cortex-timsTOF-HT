
# Import necessary libraries
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

# Define the sizes of each group and their intersection for the protein data
human_proteins = 5806
nHP_proteins = 5472
both_proteins = 5142

# Create the Venn diagram for proteins
venn = venn2(subsets = (human_proteins - both_proteins, nHP_proteins - both_proteins, both_proteins),
      set_labels = ('Human Proteins', 'nHP Proteins'),
      set_colors=('skyblue', 'yellowgreen'))

# Format the subset sizes with commas
for text in venn.subset_labels:
    text.set_text('{:,}'.format(int(text.get_text())))
venn.set_fontsize(18)
plt.show()


# Define the sizes of each group and their intersection for the peptide data
human_peptides = 78674
nhp_peptides = 83258
both_peptides = 58745

# Create the Venn diagram for peptides
venn = venn2(subsets = (human_peptides - both_peptides, nhp_peptides - both_peptides, both_peptides),
      set_labels = ('Human Peptides', 'nHP Peptides'),
      set_colors=('skyblue', 'yellowgreen'))

# Format the subset sizes with commas
for text in venn.subset_labels:
    text.set_text('{:,}'.format(int(text.get_text())))

plt.show()
