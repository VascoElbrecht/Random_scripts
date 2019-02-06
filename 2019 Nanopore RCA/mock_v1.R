# nanopore madness

setwd("~/Documents/University/2018 Guelph/2019 nano RCA/1 bioinformatics/Mock50_v1")

# import fastq
fastq <- readLines("~/Documents/University/2018 Guelph/2019 nano RCA/Raw_data_mock50_v1/FAK42214_477fdb09921842662b221840eee0fbc8d6530bea_3.fastq")

# convert to table
fastq <- data.frame("ID"=fastq[seq(1, length(fastq), 4)], "sequ"=fastq[seq(2, length(fastq), 4)], "Q"=fastq[seq(4, length(fastq), 4)], "length"= nchar(fastq[seq(2, length(fastq), 4)]), stringsAsFactors=F)


# remove sequences below 5000 bp
nrow(fastq)

fastq <- fastq[fastq$length>=5000,]
nrow(fastq)




#####
# insert plotting functions for dotplot here
#####

dir.create("dotplot")





# get matches to tag DB

data <- read.csv("out.txt", sep="\t", stringsAsFactors=F, header=F)

# count how oftenm tag was matched
temp <- data.frame(table(data$V1), stringsAsFactors=F)
temp <- temp[temp$Freq>=10,]

IDs <- as.character(temp$Var1)

dir.create("split")


for (i in 1:length(IDs)){

#get subset of 
matched <- data[data$V1==IDs[i],]
matched <- matched[order(matched$V7, decreasing=F),]

sequ <- fastq[fastq$ID==paste("@", IDs[i], sep=""),]

cat("", file=paste("split/", i, ".fastq", sep=""), append=F, sep="")

for (k in 1:(nrow(matched)-1)){


export <- c(paste(sequ$ID, "___", k, sep=""),
substr(sequ$sequ, matched$V7[k], matched$V7[k+1]),
"+",
substr(sequ$Q, matched$V7[k], matched$V7[k+1]))

rund <- round(nchar(export[2])/750)

if(rund<2){
cat(export, file=paste("split/", i, ".fastq", sep=""), append=T, sep="\n")
} else { # split long in middle

glumanda <- c(1, (nchar(export[2])/rund)*1:rund)
for (y in 1:rund){ # save split

export2 <- c(paste(sequ$ID, "___", k, sep=""),
substr(export[2], glumanda[y], glumanda[y+1]),
"+",
substr(export[4], glumanda[y], glumanda[y+1]))
cat(export2, file=paste("split/", i, ".fastq", sep=""), append=T, sep="\n")


} # split extra


} # if looop end

} # saving stuff end

} # loop end






