# nanopore madness

setwd("~/Documents/University/2018 Guelph/2019 nano RCA/1 bioinformatics/Mock50_v1")


fastafile <- list.files("../../Raw_data_mock50_v1", full.names=T)

readsperfile <- NULL

m <-1
for (m in 1:length(fastafile)){ # multi file processing!



# import fastq
fastq <- readLines(fastafile[m])

# convert to table
fastq <- data.frame("ID"=fastq[seq(1, length(fastq), 4)], "sequ"=fastq[seq(2, length(fastq), 4)], "Q"=fastq[seq(4, length(fastq), 4)], "length"= nchar(fastq[seq(2, length(fastq), 4)]), stringsAsFactors=F)


# remove sequences below 5000 bp
nrow(fastq)

fastq <- fastq[fastq$length>=5000,]
nrow(fastq)




#####
# insert plotting functions for dotplot here
#####
if(F){
dir.create("dotplot")

library("seqinr")

for(i in 1:nrow(fastq)){
tp <- unlist(strsplit(fastq$sequ[i], ""))

name <- paste("dotplot/", sub("(.*) runid.*", "\\1", fastq$ID[i]), ".pdf", sep="")

pdf(name, height=8, width=8)
dotPlot(tp, tp, wsize=68, wstep=30, nmatch=23)
dev.off()


nchar(fastq$sequ[i])
}

}






########
# run usearch
A <- system2("usearch", paste("-usearch_local ", fastafile[m], " -db REF_TAGS.txt -strand both -id 0.7 -blast6out out.txt -maxaccepts 0 -maxrejects 0 -alnout align.txt -mosaic -mincols 80", sep=""))


# get matches to tag DB

data <- read.csv("out.txt", sep="\t", stringsAsFactors=F, header=F)

# count how oftenm tag was matched
temp <- data.frame(table(data$V1), stringsAsFactors=F)
temp <- temp[temp$Freq>=10,]

nrow(temp)

IDs <- as.character(temp$Var1)

readsperfile[m] <- length(IDs)

dir.create("split")


for (i in 1:length(IDs)){

#get subset of 
matched <- data[data$V1==IDs[i],]
matched <- matched[order(matched$V7, decreasing=F),]

sequ <- fastq[fastq$ID==paste("@", IDs[i], sep=""),]

#overwrite old files
cat("", file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fastq", sep=""), append=F, sep="")
cat("", file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fasta", sep=""), append=F, sep="")


for (k in 1:(nrow(matched)-1)){


export <- c(paste(sub("(.*) runid.*", "\\1", sequ$ID), "___", k, sep=""),
substr(sequ$sequ, matched$V7[k]+58, matched$V7[k+1]+58),
"+",
substr(sequ$Q, matched$V7[k]+58, matched$V7[k+1]+58))

rund <- round(nchar(export[2])/750)

if(rund<2){
cat(export, file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fastq", sep=""), append=T, sep="\n")
#save fasta
cat(c(sub("@", ">", export[1]), export[2]), file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fasta", sep=""), append=T, sep="\n")

} else { # split long in middle

glumanda <- c(1, (nchar(export[2])/rund)*1:rund)
for (y in 1:rund){ # save split

export2 <- c(paste(export[1], "__", y, sep=""),
substr(export[2], glumanda[y]+58, glumanda[y+1]+58),
"+",
substr(export[4], glumanda[y]+58, glumanda[y+1]+58))
cat(export2, file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fastq", sep=""), append=T, sep="\n")


#save fasta
cat(c(sub("@", ">", export2[1]), export2[2]), file=paste("split/", sub("(.*) runid.*", "\\1", IDs[i]), ".fasta", sep=""), append=T, sep="\n")


} # split extra


} # if looop end

} # saving stuff end

} # loop end


########
# Make sequence alignments
# MAFFT

dir.create("alignments")

split <- list.files("split", full.names=T, pattern= ".fasta")


for(i in 1:length(split)){
A <- system2("mafft", paste("--auto ", split[i], " > ", sub("split(.*).fasta", "alignments\\1_aln.fasta", split[i]), " --thread 12", sep=""))
}


print(m)
} # end multifiles processing



write.csv(data.frame(fastafile, readsperfile), file="reads_per_4k.csv")




# genious cons




# JAMP


library("JAMP")


Empty_folder()



Cutadapt(forward="GGTCAACAAATCATAAAGAYATYGG", reverse="TAAACTTCAGGGTGACCAAARAAYCA", anchoring=F, fastq=F, bothsides=F)





