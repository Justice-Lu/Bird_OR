### Have ORF Info from Step 3d
### Have Starting M infrom from step 5d

args = commandArgs(trailingOnly=TRUE)

# Set arguments to variables 
SPP_PATH=args[1]
spp=tail(unlist(strsplit(SPP_PATH, '/')), 1) # Grabs the last directory as spp, should be Specie
setwd(SPP_PATH)


### Have coordinates of hit from name 
bed <- read.table(file.path(SPP_PATH,paste0(spp, "_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed")))
colnames(bed) <- c("Scaffold",
                    "BedStartMinus300",
                    "BedStopPlus300",
                    "bedID","evalue","strand") # these are zero based!


#### These are the results of picking ORFs in step 3
orfCoords <- read.table(file.path(SPP_PATH, "ORF_longThan_810_bp_step_2_results.olfacUniqueBlastHits.stranded.fasta.ORFStStopTable"),header=F)
colnames(orfCoords) <- c("bedID",
                         "frame",
                         "orfStart",
                         "orfStop") # these are 1 based ,"finalBedStart","finalBedEnd"

# Takes the first half of the bedID 
orfCoords$bedID <- sub("::.*", "", orfCoords$bedID)
orfCoords_0Based <- orfCoords
orfCoords_0Based$orfStart <- orfCoords$orfStart - 1

# Read the content of the file and remove 'Human_OR2J3\t'
# TODO artifact from previous file outputs consider removing from previous script. 
MCoords <- readLines(file.path(SPP_PATH, "step_5_result.pickedM.CoordinatesOfBestStart.txt"))
MCoords <- gsub("Human_OR2J3\t", "", MCoords)

MCoords <- as.data.frame(do.call(rbind, strsplit(MCoords, "\t")))
colnames(MCoords) <- c("bedID",
                       "MCoord")

MCoords$MCoord_nt <- as.integer(MCoords$MCoord) *3

# Takes the first half of the bedID 
MCoords$bedID <- sub("::.*", "", MCoords$bedID)

# Merge
bedOrf <- (merge(bed,orfCoords_0Based,by="bedID")) 
bedOrf$AACount <- (bedOrf$orfStop - bedOrf$orfStart)/3
bedOrf$bedIDPlusAA <- paste(bedOrf$bedID,"-",bedOrf$AACount,"_aa",sep="")
bedOrfM <- merge(bedOrf, MCoords,by="bedID") # ends up being dimensions of M coords
bedOrfM$orfStartCorrectM <- bedOrfM$orfStart + bedOrfM$MCoord_nt

bedOrfM$finalBedStart <- NA
bedOrfM$finalBedEnd <- NA


# for positive strand: do bedStart + orfStartCorrectM to get final start position
# and bed start + orfEnd for final end position 
# for negative strand: do bedEnd - orfSt and bedend - orfEnd.
bedOrfM[bedOrfM$strand=="+",]$finalBedStart <- bedOrfM[bedOrfM$strand=="+",]$BedStartMinus300 + bedOrfM[bedOrfM$strand=="+",]$orfStartCorrectM
bedOrfM[bedOrfM$strand=="+",]$finalBedEnd <- bedOrfM[bedOrfM$strand=="+",]$BedStartMinus300 + bedOrfM[bedOrfM$strand=="+",]$orfStop # I think this doesn't change based on M because it's all relative to the original sequence pulled out of the genome (taht is between bedstart and bedend)
bedOrfM[bedOrfM$strand=="-",]$finalBedEnd <- bedOrfM[bedOrfM$strand=="-",]$BedStopPlus300 - bedOrfM[bedOrfM$strand=="-",]$orfStart
bedOrfM[bedOrfM$strand=="-",]$finalBedStart <- bedOrfM[bedOrfM$strand=="-",]$BedStopPlus300 - bedOrfM[bedOrfM$strand=="-",]$orfStop

finalBed <- bedOrfM[,c("Scaffold","finalBedStart","finalBedEnd","bedIDPlusAA","evalue","strand")]

# Output file
write.table(finalBed, file.path(SPP_PATH, paste0(spp, "_finalORCoordinates.all.0based.bed")),quote=F,sep="\t",row.names = F,col.names = F)

cat(paste0("==== finalBed file generated ===\n",file.path(SPP_PATH, paste0(spp, "_finalORCoordinates.all.0based.bed\n\n"))))
