#integración y normalización de datos

setwd("~/Desktop/TFM/bbdd")

library(dplyr)
library(readr)
library(tidyr)
library(rvest)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(stringr)




# descargar mirbase y abrir csv

url <- "https://www.mirbase.org/results/?organism=hsa&chromosome=&start=&end="
webpage <- read_html(url)
table <- webpage %>%
  html_table(fill = TRUE) %>%
  .[[1]]  
write.csv(table, "miRBase.csv", row.names = FALSE)

datos_refseq <- read_csv("datos_refseq")
datos_ensembl <- read_csv("datos_ensembl")
datos_miRBase <- read_csv("datos_miRBase")
datos_lncipedia <- read_csv("datos_lncipedia")
datos_mircancer <- read_csv("datos_mircancer")
omim <- read_csv("omim.csv")





# prefijos mirbase

obtain_name <- function(data, column_name) {
  if (!column_name %in% colnames(data)) {
    stop(paste("La columna", column_name, "no existe en el data frame."))
  }
  
  data %>%
    mutate(prefix = substr(!!sym(column_name), 1, 7)) %>%
    filter(prefix %in% c("hsa-let", "hsa-mir")) %>%
    mutate(ncRNA_ID = case_when(
      str_starts(!!sym(column_name), "hsa-mir-") ~ paste0("MIR", toupper(str_replace_all(str_sub(!!sym(column_name), 9), "-([0-9a-zA-Z])", "\\1"))),
      str_starts(!!sym(column_name), "hsa-let-") ~ paste0("MIRLET", toupper(str_replace_all(str_sub(!!sym(column_name), 9), "-([0-9a-zA-Z])", "\\1"))),
      TRUE ~ !!sym(column_name)
    ))
}

datos_mircancer <- obtain_name(datos_mircancer, "mirId")
datos_miRBase <- obtain_name(datos_miRBase, "ncRNA_name")





# pasar ID a symbol

gene_ids <- gsub("\\.\\d+$", "", datos_ensembl$gene_id)  
symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
datos_ensembl$gene_symbol <- ifelse(is.na(datos_ensembl$gene_symbol), symbols, datos_ensembl$gene_symbol)
datos_ensembl <- datos_ensembl[!is.na(datos_ensembl$gene_symbol), ]






# estandarizar columnas


  #gene id

stand_gene <- function(data) {
  current_names <- names(data)
  gene_variants <- c("gene_symbol", "name", "Symbol")
  name_to_change <- current_names[current_names %in% gene_variants]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Gene_ID"
  }
  return(data)
}

datos_refseq <- stand_gene(datos_refseq)
datos_lncipedia <- stand_gene(datos_lncipedia)
datos_ensembl <- stand_gene(datos_ensembl)

datos_mircancer <- datos_mircancer %>%
  mutate(Gene_ID = ncRNA_ID)

datos_miRBase <- datos_miRBase %>%
  mutate(Gene_ID = ncRNA_ID)



  # ncrna id

datos_ensembl <- datos_ensembl %>%
  mutate(ncRNA_ID = Gene_ID)

datos_refseq <- datos_refseq %>%
  mutate(ncRNA_ID = Gene_ID)

datos_lncipedia <- datos_lncipedia %>%
  mutate(ncRNA_ID = Gene_ID)

stand_altncrna <- function(data) {
  current_names <- names(data)
  alt_ncrna <- c("Name", "mirID", "ncRNA Name")
  name_to_change <- current_names[current_names %in% alt_ncrna]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "ncRNA_name"
  }
  return(data)
}

datos_miRBase <- stand_altncrna(datos_miRBase)
datos_mircancer <- stand_altncrna(datos_mircancer)



  # chr

stand_chr <- function(data) {
  current_names <- names(data)
  chromosome_variants <- c("chromosome", "Chromosome", "Cromosoma", "Chromosomes")
  name_to_change <- current_names[current_names %in% chromosome_variants]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Chromosome"
  }
  return(data)
}

datos_refseq <- stand_chr(datos_refseq)
datos_lncipedia <- stand_chr(datos_lncipedia)
datos_ensembl <- stand_chr(datos_ensembl)
datos_miRBase <- stand_chr(datos_miRBase)

norm_chr <- function(data) {
  data %>%
    mutate(Chromosome = ifelse(grepl("^chr", Chromosome), Chromosome, paste0("chr", Chromosome)))
}

datos_refseq <- norm_chr(datos_refseq)
datos_ensembl <- norm_chr(datos_ensembl)
datos_miRBase <- norm_chr(datos_miRBase)



  # start

stand_start <- function(data) {
  current_names <- names(data)
  start_variants <- c("start", "Inicio", "Start", "Annotation.Genomic.Range.Start")
  name_to_change <- current_names[current_names %in% start_variants]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Start"
  }
  return(data)
}

datos_refseq <- stand_start(datos_refseq)
datos_lncipedia <- stand_start(datos_lncipedia)
datos_ensembl <- stand_start(datos_ensembl)
datos_miRBase <- stand_start(datos_miRBase)



  #end

stand_end <- function(data) {
  current_names <- names(data)
  end_variants <- c("end", "Fin", "End", "Annotation.Genomic.Range.Stop")
  name_to_change <- current_names[current_names %in% end_variants]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "End"
  }
  return(data)
}

datos_refseq <- stand_end(datos_refseq)
datos_lncipedia <- stand_end(datos_lncipedia)
datos_ensembl <- stand_end(datos_ensembl)
datos_miRBase <- stand_end(datos_miRBase)



  #estandarizar forward y reverse strand

norm_strand <- function(data) {
  current_names <- names(data)
  strand <- c("Orientation", "strand")
  name_to_change <- current_names[current_names %in% strand]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Strand"
  }
  
  data <- data %>%
    mutate(Strand = case_when(
      Strand %in% c("plus", "+", "+1") ~ "1",
      Strand %in% c("minus", "-") ~ "-1",
      TRUE ~ as.character(Strand)  # Asegurar que los valores se conviertan a caracteres
    ))
  
  return(data)
}

datos_refseq <- norm_strand(datos_refseq)
datos_ensembl <- norm_strand(datos_ensembl)
datos_miRBase <- norm_strand(datos_miRBase)



  #estandarizar deescripción

stand_desc <- function(data) {
  current_names <- names(data)
  description <- c("description")
  name_to_change <- current_names[current_names %in% description]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Description"
  }
  return(data)
}

datos_refseq <- stand_desc(datos_refseq)
datos_lncipedia <- stand_desc(datos_lncipedia)
datos_ensembl <- stand_desc(datos_ensembl)
datos_miRBase <- stand_desc(datos_miRBase)
datos_mircancer <- stand_desc(datos_mircancer)



  #estandarizar enfermedad

stand_disease <- function(data) {
  current_names <- names(data)
  disease <- c("cancer type", "Cancer")
  name_to_change <- current_names[current_names %in% disease]
  
  if (length(name_to_change) > 0) {
    index_to_change <- which(current_names %in% name_to_change) 
    names(data)[index_to_change] <- "Disease" 
  }
  
  return(data)
}

datos_refseq <- stand_disease(datos_refseq)
datos_lncipedia <- stand_disease(datos_lncipedia)
datos_ensembl <- stand_disease(datos_ensembl)
datos_miRBase <- stand_disease(datos_miRBase)
datos_mircancer <- stand_disease(datos_mircancer)

# mapear ID de OMIM para añadirlos en nuestros csv

omim$"Class ID" <- sub("http://purl.bioontology.org/ontology/OMIM/", "", omim$"Class ID")


get_OMIM_ID <- function(disease_name, omim_data) {
  disease_name_lower <- tolower(disease_name)
  preferred_labels_lower <- tolower(omim_data$`Preferred Label`)
  
  idx <- which(preferred_labels_lower == disease_name_lower)
  if (length(idx) > 0) {
    return(omim_data$`Class ID`[idx[1]])
  } else {
    return(NA)
  }
}

datos_lncipedia$OMIM_ID <- sapply(datos_lncipedia$Disease, function(x) get_OMIM_ID(x, omim))
datos_mircancer$OMIM_ID <- sapply(datos_mircancer$Disease, function(x) get_OMIM_ID(x, omim))



  # estandarizar método 

stand_method <- function(data) {
  current_names <- names(data)
  method <- c("methods")
  name_to_change <- current_names[current_names %in% method]
  
  if (length(name_to_change) > 0) {
    index_to_change <- which(current_names %in% name_to_change) 
    names(data)[index_to_change] <- "Method" 
  }
  
  return(data)
}

datos_refseq <- stand_method(datos_refseq)
datos_lncipedia <- stand_method(datos_lncipedia)
datos_ensembl <- stand_method(datos_ensembl)
datos_miRBase <- stand_method(datos_miRBase)
datos_mircancer <- stand_method(datos_mircancer)


  # estandarizar Ensembl ID

stand_ensemblid <- function(data) {
  current_names <- names(data)
  ensembl_name <- c("Ensembl.GeneIDs", "gene_id")
  name_to_change <- current_names[current_names %in% ensembl_name]
  
  if (length(name_to_change) > 0) {
    index_to_change <- which(current_names %in% name_to_change) 
    names(data)[index_to_change] <- "Ensembl_GeneID" 
  }
  
  return(data)
}

datos_refseq <- stand_ensemblid(datos_refseq)
datos_ensembl <- stand_ensemblid(datos_ensembl)


  # estandarizar muestras 

stand_sample <- function(data) {
  current_names <- names(data)
  sample <- c("sample")
  name_to_change <- current_names[current_names %in% sample]
  
  if (length(name_to_change) > 0) {
    index_to_change <- which(current_names %in% name_to_change) 
    names(data)[index_to_change] <- "Sample" 
  }
  
  return(data)
}

datos_refseq <- stand_sample(datos_refseq)
datos_lncipedia <- stand_sample(datos_lncipedia)
datos_ensembl <- stand_sample(datos_ensembl)
datos_miRBase <- stand_sample(datos_miRBase)
datos_mircancer <- stand_sample(datos_mircancer)



  # estandarizar pubmed 

stand_pubmed <- function(data) {
  current_names <- names(data)
  pm_title <- c("title", "Pubmed Article")
  pm_id <- c("pubmed id", "")
  
  name_to_change <- current_names[current_names %in% pm_title]
  if (length(name_to_change) > 0) {
    index_to_change <- which(current_names %in% name_to_change) 
    names(data)[index_to_change] <- "Pubmed_Title" 
  }
  
  name_to_change_id <- current_names[current_names %in% pm_id]
  if (length(name_to_change_id) > 0) {
    index_to_change_id <- which(current_names %in% name_to_change_id) 
    names(data)[index_to_change_id] <- "Pubmed_ID" 
  }
  
  return(data)
}

datos_refseq <- stand_pubmed(datos_refseq)
datos_lncipedia <- stand_pubmed(datos_lncipedia)
datos_ensembl <- stand_pubmed(datos_ensembl)
datos_miRBase <- stand_pubmed(datos_miRBase)
datos_mircancer <- stand_pubmed(datos_mircancer)


  #estandarizar expresión

stand_exp <- function(data) {
  current_names <- names(data)
  expression <- c("regulated", "Profile")
  name_to_change <- current_names[current_names %in% expression]
  if (length(name_to_change) > 0) {
    names(data)[names(data) == name_to_change] <- "Expression"
  }
  
  data <- data %>%
    mutate(Expression = case_when(
      Expression %in% c("up-regulated", "up", "up-Regulated") ~ "up",
      Expression %in% c("down-regulated", "down", "down-regulation") ~ "down",
      TRUE ~ as.character(Expression) 
    ))
  
  return(data)
}

datos_mircancer <- stand_exp(datos_mircancer)
datos_lncipedia <- stand_exp(datos_lncipedia)




# filtrar datos en algunos df

datos_refseq <- datos_refseq %>%
  filter(!is.na(Start) & !is.na(End))




# sustituir espacios en ID por 'barras'-'

replace_spaces <- function(data, columns) {
  for (col in columns) {
    data <- data %>%
      mutate(!!sym(col) := gsub(" ", "-", !!sym(col)))
  }
  return(data)
}

datos_refseq <- replace_spaces(datos_refseq, c("Gene_ID", "ncRNA_ID"))
datos_ensembl <- replace_spaces(datos_ensembl, c("Gene_ID", "ncRNA_ID"))
datos_lncipedia <- replace_spaces(datos_lncipedia, c("Gene_ID", "ncRNA_ID"))
datos_miRBase  <- replace_spaces(datos_miRBase, c("Gene_ID", "ncRNA_ID"))
datos_mircancer <- replace_spaces(datos_mircancer, c("Gene_ID", "ncRNA_ID"))


# obtener Start y End de genes en lugar de transcritos del gen

# ensembl

datos_ensembl_unicos <- datos_ensembl %>% distinct(Ensembl_GeneID, .keep_all = TRUE)
datos_ensembl_unicos$Ensembl_GeneID <- sub("\\..*", "", datos_ensembl_unicos$Ensembl_GeneID)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_positions <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = datos_ensembl_unicos$Ensembl_GeneID,
  mart = ensembl
)

colnames(gene_positions) <- c("Ensembl_GeneID", "Start", "End")
gene_positions$Ensembl_GeneID <- sub("\\..*", "", gene_positions$gene_id)
gene_positions$Ensembl_GeneID <- as.character(gene_positions$Ensembl_GeneID)

datos_ensembl_final <- datos_ensembl_unicos %>%
  left_join(gene_positions, by = "Ensembl_GeneID")

datos_ensembl <- datos_ensembl_final %>%
  dplyr::select(
    transcript_id, Chromosome, Start = Start.y, End = End.y, Strand, Ensembl_GeneID, gene_biotype,
    transcript_biotype, Gene_ID, Description, ncRNA_ID
  )


#refseq

datos_refseq_unicos <- datos_refseq %>% distinct(Ensembl_GeneID, .keep_all = TRUE)
datos_refseq_unicos$Ensembl_GeneID <- sub("\\..*", "", datos_refseq_unicos$Ensembl_GeneID)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_positions <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = datos_refseq_unicos$Ensembl_GeneID,
  mart = ensembl
)

colnames(gene_positions) <- c("Ensembl_GeneID", "Start", "End")
gene_positions$Ensembl_GeneID <- sub("\\..*", "", gene_positions$Ensembl_GeneID)
gene_positions$Ensembl_GeneID <- as.character(gene_positions$Ensembl_GeneID)

datos_refseq_final <- datos_refseq_unicos %>%
  left_join(gene_positions, by = "Ensembl_GeneID")

datos_refseq <- datos_refseq_final %>%
  dplyr::select(
    NCBI.GeneID, 
    Gene_ID, 
    Description, 
    Transcripts, 
    Chromosome, 
    Ensembl_GeneID, 
    Annotation.Genomic.Range.Accession, 
    Start = Start.y, 
    End = End.y, 
    Strand, 
    Synonyms, 
    ncRNA_ID
  )


                                
# guardar los datos normalizados y estandarizados

write.csv(datos_refseq, "datos_refseq", row.names = FALSE)
write.csv(datos_ensembl, "datos_ensembl", row.names = FALSE)
write.csv(datos_miRBase, "datos_miRBase", row.names = FALSE)
write.csv(datos_lncipedia, "datos_lncipedia", row.names = FALSE)        
write.csv(datos_mircancer, "datos_mircancer", row.names = FALSE)


