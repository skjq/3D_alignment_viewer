# Name: 3Dav_main.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 10/04/2015
# Desc: main script to produce graph of multiple sequence alignment, will need a 
#       pdb file and fasta format protein MSA

source('3Dav_header.R')

## data loading steps
p = c('select the protein multiple sequence alignment file in FASTA format, it', 
      'is assumed that the first sequence in the MSA is the one with pdb sequence')
print(paste(p, collapse = ' '))
csAlFile = file.choose()

p = c('select the pdb file')
print(p)
csPdbFile = file.choose()

# load the pdb file
# create a CChain object to load the chain data
oChain = CChain(pdb.file = csPdbFile)
# if no errors in creating chain, load chain data
c = getChain(oChain)
p = getPosition(oChain)
r = getResidueNumber(oChain)
s = getSequence(oChain)

pr = c('unique chains found', unique(c), 
       'select the chain used in the alignment')
print(paste(pr, collapse = ' '))
cChain = readLines(n = 1)
# select position, residue number and sequence of selected chain
f = (c == cChain)
cvSeq.pdb = s[f]
ivRes = r[f]
mPos = p[f,]

# read the alignment file
oPwa = readAAMultipleAlignment(csAlFile)

# create mapping object of class CMapChainToFasta
# the sequence from alignment file has spaces in it 
# and some residues may not be present in the chain sequence
# x, y, z position of alpha carbon of each pdb sequence residue needs to be
# mapped to the corresponding position in the alignment sequence
cvSeq.fasta = unlist(strsplit(as.character(oPwa)[1], ''))
oMapSeq = CMapChainToFasta(cvSeq.pdb, cvSeq.fasta, mPos)

# remove the regions in the alignment matrix with no residues i.e. spaces
# get matrix alignment
mAl = as.matrix(oPwa)
colnames(mAl) = 1:ncol(mAl)
# reduce alignment matrix to only regions with coordinates i.e. no spaces
mAl = mAl[,getResidueIndex(oMapSeq)]

# now we are ready to make graphs
# repeat this function, depending on how many sequences in the alignment file
lig = vector('list', length = nrow(mAl))
names(lig) = rownames(mAl)
# make the graphs
for (i in 1:length(lig)){
  lig[[i]] = f_oIGsequenceToGraph(mAl[i,], colnames(mAl), oMapSeq)
}

# the weights of the graph are normally distributed
# and the higher the weight, the closer the residues are to each other
sapply(seq_along(lig), function(x) summary(E(lig[[x]])$weight))
plot(density(E(lig[[1]])$weight), main='distribution of weight')

p.old = par(mar=c(1,1,1,1))

for (i in 1:length(lig)){
  ig = lig[[i]]
  w = E(ig)$weight
  f = which(w < quantile(w, 0.75))
  ig.1 = delete.edges(ig, edges = f)
  n = make.names(paste(V(ig.1)$residue_name, V(ig.1)$name))
  plot(ig.1, vertex.size=0.1, vertex.label=n, layout=layout.fruchterman.reingold,
       main=names(lig)[i])
  lig[[i]] = ig.1
}

# intersect the 3 graphs to make a common graph
ig = graph.intersection(lig)
E(ig)$weight = E(ig)$weight_1
plot(ig, vertex.size=0.1, layout=layout.fruchterman.reingold)

dir.create('Results', showWarnings = F)
write.graph(ig, file='Results/ig_test.graphml', format = 'graphml')
