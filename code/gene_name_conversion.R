setwd("~/Documents/Projects/why-pathway-methods-work/code")
library(biomaRt)

OUTFILE='../data/omnipath/omnipath_uni_hgnc.csv'
data=read.csv('../data/omnipath/all_interactions.csv',
              header=T,stringsAsFactors=FALSE,sep=',',row.names = 1)
genelist=unique(c(data$target,data$source))


#possibilities: "hsapiens_gene_ensembl", ="mmusculus_gene_ensembl", "rnorvegicus_gene_ensembl"
species='hsapiens_gene_ensembl'
#possibilities: 'hgnc_symbol', "mgi_symbol","rgd_symbol"
#possibilities: 'entrezgene','ensembl_gene_id','uniprotswissprot','uniprotsptrembl'
from='uniprotswissprot'
to='hgnc_symbol'

mart = useMart("ensembl", dataset = species)

genelist=getBM(values=genelist,attributes = c(from,to), 
                filters = from,mart = mart)

write.csv(genelist,file=OUTFILE)
