# Name: 3Dav_main.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 10/04/2015
# Desc: main script to produce graph of multiple sequence alignment, will need a 
#       pdb file and fasta format protein MSA

source('3Dav_header.R')

## data loading steps
p = c('select the protein multiple sequence alignment file in FASTA format, it', 
      'is assumed that the all the sequences in the MSA have a pdb file')
print(paste(p, collapse = ' '))
csAlFile = file.choose()
# read the alignment file
oPwa = readAAMultipleAlignment(csAlFile)
len = dim(oPwa)[1]
# list to hold graphs
lig = vector('list', length = len)
names(lig) = rownames(oPwa)

# create the graphs
for (i in 1:len){
  p = c('select the pdb file')
  print(paste(p, 'for sequence', rownames(oPwa)[i]))
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
  # create mapping object of class CMapChainToFasta
  # the sequence from alignment file has spaces in it 
  # and some residues may not be present in the chain sequence
  # x, y, z position of alpha carbon of each pdb sequence residue needs to be
  # mapped to the corresponding position in the alignment sequence
  cvSeq.fasta = unlist(strsplit(as.character(oPwa)[i], ''))
  oMapSeq = CMapChainToFasta(cvSeq.pdb, cvSeq.fasta, mPos)
  # remove the regions in the alignment matrix with no residues i.e. spaces
  # get matrix alignment
  mAl = as.matrix(oPwa)
  colnames(mAl) = 1:ncol(mAl)
  # reduce alignment matrix to only regions with coordinates i.e. no spaces
  mAl = mAl[,getResidueIndex(oMapSeq)]
  # now we are ready to make graphs
  # repeat this function, depending on how many sequences in the alignment file
  # make the graphs
  lig[[i]] = f_oIGsequenceToGraph(mAl[i,], colnames(mAl), oMapSeq)
}

# the weights of the graph are normally distributed
# and the higher the weight, the closer the residues are to each other
sapply(seq_along(lig), function(x) summary(E(lig[[x]])$weight))
plot(density(E(lig[[1]])$weight), main='distribution of weight')

p.old = par(mar=c(1,1,1,1))
print('Choose cutoff distance ')
c = as.numeric(readLines(n = 1))

for (i in 1:length(lig)){
  ig = lig[[i]]
  w = E(ig)$weight
  f = which(w < (max(w)-c))
  ig.1 = delete.edges(ig, edges = f)
  n = make.names(paste(V(ig.1)$residue_name, V(ig.1)$name))
  plot(ig.1, vertex.size=0.1, vertex.label=n, layout=layout.fruchterman.reingold,
       main=names(lig)[i])
  lig[[i]] = ig.1
}

# intersect the graphs to make a common graph
ig = graph.intersection(lig)
E(ig)$weight = E(ig)$weight_1
plot(ig, vertex.size=0.1, layout=layout.fruchterman.reingold)

dir.create('Results', showWarnings = F)
write.graph(ig, file='Results/ig_test.graphml', format = 'graphml')

# further convert the graphs using bipartite graph weighting strategy
# using CGraph class
# get the adjacency matrix
mInc = get.adjacency(ig, 'both', sparse=FALSE)
# set diagonals to 1 as to have connection with itself
diag(mInc) = 1
# change rownames as those will be vertex of type 2 i.e. type FALSE
rownames(mInc) = paste(rownames(mInc), rownames(mInc), sep='.')
# create bipartite graph
ig.p = graph.incidence(mInc)
# sanity check
if (!is.bipartite(ig.p)) stop('Graph is not bipartite, check error')
# create CGraph object to assign weights based on interestingness
ob = CGraph(oGraph = ig.p)
ig.p = getProjectedGraph(ob)

# set new weight values
E(ig.p)$neighbours = E(ig.p)$weight
E(ig.p)$weight = E(ig.p)$ob_to_ex
w = E(ig.p)$weight
summary(w)
quantile(w)
par(p.old)
hist(w)
# choose a cutoff by modelling the distribution shape
# it appears that the distribution follows a power law?
# taking square root means we can fit a poisson distribution
w2 = sqrt(w)
r = round(range(w2))
s = seq(r[1]-0.5, r[2]+0.5, by = 1)
hist(w2, prob=T, breaks=s, main='distribution of obs to exp ratios', 
     xlab='square root obs to exp ratio', ylab='')
dp = dpois(r[1]:r[2], lambda = mean(w2))
dn = dnbinom(r[1]:r[2], size = mean(w2), mu = mean(w2))
lines(r[1]:r[2], dp, col='red')
lines(r[1]:r[2], dn, col='blue')
legend('topright', legend = c('poi', 'nbin'), fill = c('red', 'blue'))

print(paste('Choose cutoff point looking at histogram, values less than cutoff will',
      'be removed'))
c = as.numeric(readLines(n=1))

f = which(w2 < c)
ig.p = delete.edges(ig.p, f)
d = degree(ig.p)
ig.p = delete.vertices(ig.p, which(d == 0))
par(mar=c(1,1,1,1))
plot(ig.p, vertex.size=0.1, layout=layout.fruchterman.reingold)

# save the graph for examining using cytoscape
write.graph(ig.p, file='Results/ig_test_bp.graphml', format = 'graphml')

# intersect the 2 graphs to make final graph
ig.f = graph.intersection(ig, ig.p)
# set weight parameters
E(ig.f)$distance = E(ig.f)$weight_1
E(ig.f)$weight = E(ig.f)$ob_to_ex

# save the graph for examining using cytoscape
write.graph(ig.f, file='Results/ig_final.graphml', format = 'graphml')



