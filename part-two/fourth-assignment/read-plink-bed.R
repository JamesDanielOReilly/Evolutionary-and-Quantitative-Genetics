read.plink = function(root) { # function to read in PLINK .bed files, giving the genotypes as 0/1/2/NA, each row is one individual -- beware the script will try to read in the whole data in the bed file, so make sure it is small enough to be read in in one go
  bed.file = paste(root, '.bed', sep = '')   
  bed.file.size = file.info(bed.file)$size # the bed file size in bytes
  sample.size = dim(read.table(paste(root, '.fam', sep = '')))[1] # number of individuals
  snp.size = ceiling(sample.size/4) # no. bytes storing the data for 1 SNP for all individuals; each byte stores 4 people
  n.snps = round((bed.file.size-3)/snp.size) # the total .bed file size is 3 plus size of ids x snps combo hence removing 3
  bin.connection = file(bed.file, 'rb') # opens a connection with .bed file
  test.bytes = readBin(bin.connection, what = "raw", n = 3) # bytes 1 and 2 store infor about file type, byte 3 about array type
  if(!identical(as.character(test.bytes), c('6c', '1b', '01'))) {
    stop('BED file not a v0.99 SNP-major BED file, please re-encode the data as v0.99 SNP-major file')
  }
  genotypes = matrix(ncol = n.snps, nrow = sample.size) # to store the genos (obviously)
  for(i in 1:n.snps) {
    r.bin.snp = readBin(bin.connection, what = 'raw', n = snp.size)
    bin.snp = matrix(as.numeric(rawToBits(r.bin.snp)), ncol = 2, byrow = TRUE)[1:sample.size,]
    genotypes[,i] = bin.snp[,1] + bin.snp[,2] - 10 * ((bin.snp[,1] == 1) & (bin.snp[,2] == 0))
  }
  genotypes[genotypes == -9] = NA
  snp.names = read.table(paste(root, '.bim', sep = ''))[,2] # the bim file has the snps names in the 2nd column
  colnames(genotypes) = snp.names # give the names to the snps
  close(bin.connection) # needed to avoid pointless warning
  return(genotypes) # must be after closing the connection or the pointless warning reappears
}

