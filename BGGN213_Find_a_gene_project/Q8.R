library(bio3d)
library(gplots)
# Read fasta multiple alignment file
fasta <- read.fasta("muscle_alignment2.fst")

# Make an identity matrix and heatmap it
identity <- seqidentity(fasta)

heatmap(identity, margins = c(11,5), symm = T, col = heat.colors(256))
dim(identity)

# Find a consensus sequence from the multiple alignement and blast it
con <- consensus(fasta)
seq <- paste(con$seq, collapse = '')
blast <- blast.pdb(seq)
plot.blast(blast)

# Pick out the top three unique hits
ids <- blast$hit.tbl$subjectids[c(1,2,4)]
#blast$hit.tbl$identity

# Annotate more information
pdb.annotate("4JVW", anno.terms = c("experimentalTechnique", "resolution", "source"))
pdb.annotate("1R70", anno.terms = c("experimentalTechnique", "resolution", "source"))
pdb.annotate("1IGA", anno.terms = c("experimentalTechnique", "resolution", "source"))

# Example sequence
ex <- "TKFAKLSCLVTNLATYDTLNISWSSKSGEPLETNTKIMESHPNGTFSAVGVASVCMEDWDNRKEFVCTVTHRDLPSPQKKFISKPNEVAKHPPAVYLLPPAREQLILRESATVTCLVKGFSPADIFVQWLQRGQPLSSDKYVTSAPMPEPGAPGLYFTHSILTVTEEEWNSGETYTCVVGHEALPHMVTERTVDKSTGKPTLYNVSLIMSDTGGTCYPCRSTRQALGVELLCVCKLTMSEDVAFYKKKKKKX"
#blast.pdb(ex)

