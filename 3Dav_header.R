# 3Dav_header.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 06/04/2015
# Desc: header file for the project 3D_alignment_viewer loading libraries, setting global
#       variables and global classes and functions


library(methods)


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

