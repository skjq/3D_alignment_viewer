# Name: 3Dav_header.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 06/04/2015
# Desc: header file for the project 3D_alignment_viewer loading libraries, setting global
#       variables and global classes and functions


library(methods)
library(igraph)
library(Biostrings)


#### Classes and functions
###### class CChain
# Desc: class to hold sequences, chains, residue ids and 3d position of each residue
#       from calculated via each CA ATOM line of pdb file and has generic accessor
#       functions
# Usage: create by calling CChain constructor and PDB file name/path
setClass('CChain', slots=list(dfAtoms='data.frame'))

# constructor
CChain = function(pdb.file){
  ## private functions
  # Name: f_cvGetCAatomFromPDBFile
  # Args: pdb file name
  # Rets: character string vector of ATOM Lines of pdb fle 
  # Desc: takes a pdb file name, returns CA ATOM lines
  f_cvGetCAatomFromPDBFile = function(pdb.file){
    # open file
    infile = file(pdb.file, 'rt')
    # read 1000 lines at a time
    input=readLines(infile, n=-1)
    # check for CA atom line
    i = grep('^ATOM\\s+\\d+\\s+CA{1}.+', input)
    # error checking
    if (length(i) ==0) stop(paste('No ATOM Lines found in file', pdb.file))
    # close file
    close(infile)
    return(input[i])
  } # function
  
  # Name: f_csRemoveWhiteSpace
  # Args: character string with possible white space
  # Rets: string with white space removed from beginning and end
  # Desc: removes white space from beginning and end of string
  f_csRemoveWhiteSpace <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  # Name: f_lParseAtomRecordFromPDB
  # Args: character string with atom line from pdb file
  # Rets: list with each record from atom line split
  # Desc: takes an ATOM record line from the pdb file, extracts each record from the
  #       columns of the strings according to 
  #       http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
  #       and assigns it to a list and returns the list
  f_lParseAtomRecordFromPDB = function(st)
  {
    st = unlist(strsplit(st, '', perl=T))  
    lRet = vector('list')
    lRet$Atom = f_csRemoveWhiteSpace(paste(na.omit(st[1:4]), collapse = ''))
    lRet$AtomSrNo = as.integer(paste(na.omit(st[7:11]), collapse = ''))
    lRet$AtomName = f_csRemoveWhiteSpace(paste(na.omit(st[13:16]), collapse = ''))
    lRet$AltLoc = f_csRemoveWhiteSpace(paste(na.omit(st[17]), collapse = ''))
    lRet$ResName = f_csRemoveWhiteSpace(paste(na.omit(st[18:20]), collapse = ''))
    lRet$Chain = f_csRemoveWhiteSpace(paste(na.omit(st[22]), collapse = ''))
    lRet$ResSeqNo = as.integer(paste(na.omit(st[23:26]), collapse = ''))
    lRet$InsCode = f_csRemoveWhiteSpace(paste(na.omit(st[27]), collapse = ''))
    lRet$X = as.numeric(paste(na.omit(st[31:38]), collapse = ''))
    lRet$Y = as.numeric(paste(na.omit(st[39:46]), collapse = ''))
    lRet$Z = as.numeric(paste(na.omit(st[47:54]), collapse = ''))
    lRet$Occupancy = as.numeric(paste(na.omit(st[55:60]), collapse = ''))
    lRet$Temperature = as.numeric(paste(na.omit(st[61:66]), collapse = ''))
    lRet$Temperature = as.numeric(paste(na.omit(st[61:66]), collapse = ''))
    lRet$SegmentID = f_csRemoveWhiteSpace(paste(na.omit(st[73:76]), collapse = ''))
    lRet$Element = f_csRemoveWhiteSpace(paste(na.omit(st[77:78]), collapse = ''))
    lRet$Charge = f_csRemoveWhiteSpace(paste(na.omit(st[79:80]), collapse = ''))
    return(lRet)
  }
  
  
  ## end private functions
  # read CA atoms from pdb file
  cvAtoms = f_cvGetCAatomFromPDBFile(pdb.file)
  dfAtoms = t(sapply(seq_along(cvAtoms), function(x) f_lParseAtomRecordFromPDB(cvAtoms[x])))
  dfAtoms = data.frame(dfAtoms)
  # some residues may have duplicate positions, so use the first one
  # this can be checked for each chain
  c = unique(as.character(dfAtoms$Chain))
  dfAtoms.2 = NULL
  for (i in 1:length(c)){
    # check the residues for each chain
    # use selected chain to subset the dataframe
    df = dfAtoms[as.character(dfAtoms$Chain) %in% c[i], ]
    # check if the residues in this chain are duplicated
    f = !duplicated(as.numeric(df$ResSeqNo))
    df = df[f,]
    dfAtoms.2 = rbind(dfAtoms.2, df)
  }
  new('CChain', dfAtoms=dfAtoms.2)
} # end constructor

# accessor functions
setGeneric('getChain', function(obj)standardGeneric('getChain'))
setMethod('getChain', signature = 'CChain', definition = function(obj){
  ret = obj@dfAtoms
  return(as.character(ret$Chain))
})

setGeneric('getPosition', function(obj)standardGeneric('getPosition'))
setMethod('getPosition', signature = 'CChain', definition = function(obj){
  ret = obj@dfAtoms
  m = cbind(as.numeric(ret$X), as.numeric(ret$Y), as.numeric(ret$Z))
  return(m)
})

setGeneric('getSequence', function(obj)standardGeneric('getSequence'))
setMethod('getSequence', signature = 'CChain', definition = function(obj){
  # check if biostrings library present
  if (!require(Biostrings)) stop('Biostrings Library required')
  cvRes = as.character(obj@dfAtoms$ResName)
  cvRes.2 = sapply(seq_along(cvRes), function(x){
    # get the single letter residue codes
    i = grep(cvRes[x], AMINO_ACID_CODE, ignore.case = T)
    names(AMINO_ACID_CODE[i])
    # NOTE: Add non-standard residue check later
  })
  return(cvRes.2)
})


setGeneric('getResidueNumber', function(obj)standardGeneric('getResidueNumber'))
setMethod('getResidueNumber', signature = 'CChain', definition = function(obj){
  return(as.numeric(obj@dfAtoms$ResSeqNo))
})

############ end class CChain

######### Class CMapChainToFasta
# Name: CMapChainToFasta
# Desc: the sequence from a pdb file along with its coordinates need to be mapped 
#       to the sequence from the fasta alignment file. when the fasta alignment sequence
#       is taken from an alignment file, there are usually spaces in the forms of - in 
#       the sequences, which means that mapping the sequence from the pdb file (which will
#       not have any spaces) to the same sequence from the fasta alignment file needs to be
#       done
# Args: requires an array of characters for the sequence from the pdb file i.e. arChain.
#       arFasta.seq: same sequence but from fasta alignment file with spaces
#       iCoordinates: coordinates for the arChain sequence from pdb file
# Rets: an object of the class that contains the mapped data
setClass('CMapChainToFasta', slots=list(map='array', pos='matrix'))
# constructor
CMapChainToFasta = function(arChain, arFasta.seq, mCoordinates) {
  # align the 2 sequences using biostrings class
  if (!require(Biostrings)) stop('Bioconductor library Biostrings required')
  # align pdb sequence to fasta sequence with spaces i.e. dashes -
  pwa = pairwiseAlignment(paste(arChain, collapse = ''), 
                          paste(arFasta.seq, collapse = ''))
  # create a matrix of alignment object with chain seq and fasta seq
  m = as.matrix(pwa)
  names(arFasta.seq) = NULL
  m = rbind(m, arFasta.seq)
  # find the positions of the spaces i.e. - in the alignment
  i = grep('-', m[1,])
  x = 1:ncol(m)
  # x has all the coordinates of the columns of matrix
  # i has the columns that have a -
  # subtracting i from x will give columns that have a sequence
  s.len = start(pattern(pwa)):end(pattern(pwa))
  # if perfect match with no spaces then grep returned index will be of length 0
  # in that case no missing - unmapped data, aligns perfectly
  if (length(i) > 0) x.map = x[-i] else x.map = x
  # iCoordinates maps to rows of position matrix i.e. mCoordinates
  iCoordinates = s.len
  ar = matrix(c(s.len, x.map, iCoordinates), nrow = 3, byrow = T, 
              dimnames = list(c('seq', 'map', 'coord'), NULL))
  # some residues from the pdb sequence may be aligned to gaps
  # and those will also need to be removed from this matrix as those
  # residues do not exist in the fasta sequence 
  x = ar['map',]
  # check if any of these positions in the fasta sequence
  # have a - which means there is no residue there and hence no coordinate
  i = grep('-', arFasta.seq[x])
  ar.in = 1:ncol(ar)
  # remove these indices from the mapping matrix
  if (length(i) > 0) ar.in = ar.in[-i]
  ar = ar[,ar.in]
  new('CMapChainToFasta', map=ar, pos=mCoordinates[ar.in,])  
} # end constructor

## accessor functions
# Name: getResidueIndex
# Args: object of class CMapChainToFasta
# Rets: gives the index numbers for residues that have a coordinate
setGeneric('getResidueIndex', function(obj) standardGeneric('getResidueIndex'))
setMethod('getResidueIndex', signature = 'CMapChainToFasta', definition = function(obj){
  m = obj@map['map',]
  return(round(m,0))
})

# Name: getPosition
# Args: object of class CMapChainToFasta
# Rets: gives the positions for these residues i.e. x, y, z coordinates
## this genereic already exists for CChain class
#setGeneric('getPosition', function(obj) standardGeneric('getPosition'))
setMethod('getPosition', signature = 'CMapChainToFasta', definition = function(obj){
  m = obj@pos
  return(round(m,2))
})

######### End Class CMapChainToFasta

######### Functions
# Name: f_oIGsequenceToGraph
# Desc: Converts the amino acid sequence into a graph object with each edge weighted
#       depending on how close the 2 corresponding vertices are to each other
# Args: arSeq - an array of characters of amino acid residues
#       arVertexNames - an array of residue ids that will be assigned to the vertex names
#       this should be index numbers of residues in alignment matrix after mapping of
#       pdb and fasta sequences using the class CMapChainToFasta
#       oMap - object of CMapChainToFasta class created earlier
# Rets: an object of class igraph
f_oIGsequenceToGraph = function(arSeq, arVertexNames, oMap){
  if (!require(igraph)) stop('igraph library required')
  # sequence - array of characters
  s = arSeq
  # array of vertex names - these will be numbers which correspond
  # to the positions of the residues in the alignment matrix
  # that are left behind after doing the mapping step
  r = arVertexNames
  # create the diatance matrix OR adjacency matrix - square matrix
  m = as.matrix(dist(getPosition(oMap)))
  rownames(m) = r
  colnames(m) = r
  # invert the weights, things that are farthest will have the smallest number
  # others that are closest to each other will have the most remaining behind
  # e.g. 10 - 5 = 5 no change in weight
  # 10 - 8 = 2 as these 2 are far away, they will have a weight of 2 remaining etc
  m = max(m) - m
  # no connections or weights on diagonals
  diag(m) = 0  
  ig = graph.adjacency(m, 'min', weighted = T)
  V(ig)$residue_name = s
  return(ig)
}

##### Class CGraph
# Name: Class CGgraph
# Desc: assigns weights to one mode projection of graphs based on observed to expected probabilities of 
#       vertices of the first kind i.e. with value TRUE using the igraph library
#       Zweig, K. A., & Kaufmann, M. (2011). A systematic approach to the one-mode projection of 
#       bipartite graphs. Social Network Analysis and Mining (Vol. 1, pp. 187–218). 
#       doi:10.1007/s13278-011-0021-0

# declaration
setClass('CGraph', slots=list(ig='ANY', r='numeric', f='logical', ig.p='ANY'))

# object constructor
CGraph = function(oGraph){
  # check if igraph library present
  if (!require(igraph)) stop('R library igraph required')
  # check if graph is bipartite
  if (!is.bipartite(oGraph)) stop('Graph is not bipartite')
  #### internal private functions
  # processing steps - called by constructor  
  # assign probabilities to vertex of first kind
  # Name: CGraph.assign.marginal.probabilities
  # Desc: assigns probabilities to each vertex of the first kind (TRUE) 
  #       based on how many times it is connected to the vertex of the 
  #       second kind i.e. degree(V1) / (total number of V-type2)
  # Args: internal function - object of CGraph class
  CGraph.assign.marginal.probabilities = function(obj){
    # vertex of the first kind will be assigned probabilities
    # based on their relations with the vertices of the second kind
    # flag to identify vertex types
    f = V(obj@ig)$type
    d = degree(obj@ig)
    d = d[f]
    # r is the total numbers of vertices of the second kind
    r = sum(!f)
    p = d/r
    V(obj@ig)[f]$prob_marginal = p
    obj@r = r
    obj@f = f
    return(obj)
  }
  
  # Name: CGraph.project
  # Desc: assigns a level of interestingness/leverage or observed to expected ratio to 
  #       each edge after graph projection on the vertex of first kind i.e. type = TRUE 
  #       Observed frequency = weight of edge / (total number of vertices of second type)
  #       i.e. how many shared vertices of type 2 are between the 2 type 1 vertices
  #       Expected frequency = how many times we expect to see them based on their 
  #       joint probability under assumption of independence. 
  #       (marginal.prob of V1 * marginal.prob of V2)
  # Args: called internally no need to do it externally, 
  #       will project on vertex with TYPE=TRUE
  CGraph.project = function(obj){
    # project the graph in one dimension and
    # assign weights based on observed to expected ratios
    g.p = bipartite.projection(obj@ig, which = 'TRUE')
    # get the matrix with rows representing each edge
    m = get.edgelist(g.p)
    w = E(g.p)$weight
    # calculate observed ratio
    # weight / r
    ob = w / obj@r
    # calculate expected 
    mExp = cbind(V(g.p)[m[,1]]$prob_marginal, V(g.p)[m[,2]]$prob_marginal)
    ex = mExp[,1] * mExp[,2]
    E(g.p)$observed = ob
    E(g.p)$expected = ex
    E(g.p)$ob_to_ex = ob / ex
    obj@ig.p = g.p
    return(obj)
  }
  
  ####
  # create the object
  g = new('CGraph', ig=oGraph, r = 0, f= F, ig.p=NULL)
  # assign marginal probabilities
  g = CGraph.assign.marginal.probabilities(g)
  # assign weights on one mode projection
  g = CGraph.project(g)
  return(g)
}


# data acccessor functions
setGeneric('getBipartiteGraph', function(obj)standardGeneric('getBipartiteGraph'))
setMethod('getBipartiteGraph', signature = 'CGraph', definition = function(obj){
  return(obj@ig)
})

setGeneric('getProjectedGraph', function(obj)standardGeneric('getProjectedGraph'))
setMethod('getProjectedGraph', signature = 'CGraph', definition = function(obj){
  return(obj@ig.p)
})




