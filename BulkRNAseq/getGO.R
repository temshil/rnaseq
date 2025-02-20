#Get GO
getGO <- function(org.db, keys, keytype) {

require(GO.db)

temp.GO <- select(org.db, keys = keys, keytype=keytype, columns = c("GO","ONTOLOGY"))
temp.GO <- temp.GO[!is.na(temp.GO$GO),]
temp.GO <- cbind(temp.GO, select(GO.db, key=temp.GO$GO,keytype="GOID",columns="TERM"))
temp.GO$idterm <- paste(temp.GO$GO, temp.GO$TERM, sep=" ")
temp.GO2 <- tapply(temp.GO$idterm, list(temp.GO[,1], temp.GO$ONTOLOGY), paste, collapse="; ")
return(data.frame(temp.GO2))
}
