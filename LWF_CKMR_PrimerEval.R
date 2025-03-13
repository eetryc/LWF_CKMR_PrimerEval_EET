library(vcfR)
library(tidyverse)
library(adegenet)
#Load primers
primers_nt = read.csv("primers_notail.csv")
head(primers_nt)

primers_nt = primers_nt %>% 
  dplyr::select(!`Sequence`) %>% 
  distinct()



#Load catalog windows and radtag info
catalog_winds = readRDS("final_catalog_windows_seqs_2filter.RDS")
head(catalog_winds)

load("radtag_positions.rdata")
head(radtag_positions)

#Filter catalog info by primer
catalog_winds.filt = catalog_winds %>% filter(primer3_name %in% primers_nt$Name)

#Append to the primer info
primers.radtags = primers_nt %>% 
  left_join(catalog_winds, by = c("Name" = "primer3_name"))


#Filter the radtag positions to only those radtags
radtag_positions.filt = radtag_positions %>% filter(RADtag %in% primers.radtags$RADtag)

#Load the vcf file
snps.vcf = read.vcfR("populations_filtered_miss065.snps.vcf")
head(snps.vcf)
range(catalog_winds$Flanking1_Start)

#Create a new file for the filtered VCF (start by copying over the full file so all the important bits are there)
snps.vcf.filt = snps.vcf

#Then filter the original to only contain IDs in the reduced set(RADtag:position:strand)
snps.vcf.filt@fix = snps.vcf@fix[which(as.data.frame(snps.vcf@fix)$ID %in% radtag_positions.filt$ID),]
snps.vcf.filt@gt = snps.vcf@gt[which(as.data.frame(snps.vcf@fix)$ID %in% radtag_positions.filt$ID),]
#Inspect the result
snps.vcf.filt




snps.genind.filt <- vcfR2genind(snps.vcf.filt)




#use key to add matching sample names
key=read.csv("LWF_key.csv")

modified_names <- indNames(snps.genind.filt) %>%
  as_tibble() %>%
  separate(value, into = c("indiv", NA), sep = "_") %>%
  pull(indiv)  # Convert back to a vector

#Match the modified names to key$Run
matching_indices <- match(modified_names, key$Run)

# Replace only matched names
new_names <- modified_names  # Start with modified names
new_names[!is.na(matching_indices)] <- key$Sample.Name[matching_indices[!is.na(matching_indices)]]

#  Assign correct names back to the genind object
indNames(snps.genind.filt) <- new_names

# See if it worked...
indNames(snps.genind.filt) %>% 
  as_tibble() %>% 
  print(n =5000) 



#read in pop data to match
githubinstall::jaredhomola/MCGLfxns
population=read.csv("lwfRubiasBaseline_10242023.csv") %>% 
  as_tibble() %>% 
  dplyr::select(collection,indiv)

pop.genind=snps.genind.filt[pop$indiv,]


pop_data <- data.frame(indiv = indNames(snps.genind.filt), population = pop(snps.genind.filt))
pop.genind <- snps.genind.filt[pop_data$indiv, ]

# Extract only the individuals that exist in pop_data$indiv
matching_inds <- indNames(snps.genind.filt) %in% pop_data$indiv
pop.genind <- snps.genind.filt[matching_inds]

#match the ind to populations
population %>% 
  filter(indiv %in% rownames(pop.genind@tab)) %>% 
  arrange(match(indiv, rownames(pop.genind@tab)))


pop_vector <- population$collection[matching_indices]

#Pop as factor change
pop(snps.genind.filt)<-as.factor(pop_vector)

snps.genind.filt.n <- snps.genind.filt[!is.na(snps.genind.filt$pop), ]

#PCA
mcglPCA(snps.genind.filt)
#error: "Error in select(., pop) : unused argument (pop)"





#Added dplyr:: to select functions
mcglPCA2= function(dat) {
  x <- tab(dat,
           freq = TRUE,
           NA.method = "mean")
  
  pca <- dudi.pca(
    x,
    center = TRUE,
    scale = FALSE,
    nf = 4,
    scannf = FALSE
  )
  
  popNames <- dat@pop %>%
    as_tibble() %>%
    rename(pop = value)
  
  pca.df <- pca$li %>%
    as_tibble() %>%
    mutate(pop = popNames$pop)
  
  
  centroids <- pca.df %>%
    dplyr::select(pop) %>%
    right_join(aggregate(cbind(Axis1, Axis2) ~ pop,
                         pca.df,
                         mean)) %>%
    distinct(.keep_all = TRUE)
  
  ggplot(pca.df,
         aes(x = Axis1,
             y = Axis2,
             color = pop)) +
    geom_point(alpha = 0.4) +
    stat_ellipse() +
    geom_point(data = centroids,
               size = 6,
               alpha = 0.8) +
    geom_text_repel(
      data = centroids,
      fontface = "bold",
      aes(label = pop),
      size = 5,
      force = 5,
      force_pull = 0.1,
      max.overlaps = 100
    ) +
    ylab("Axis 2") +
    xlab("Axis 1") +
    theme_classic(base_size = 18) +
    theme(legend.position = "none") +
    NULL
}

mcglPCA2(snps.genind.filt)


#to test only a few 
library(poppr)

subset.genind <- snps.genind.filt %>% popsub(sublist = c("Fox River L1", "Muskegon", "Ingalls Point and Ingalls Bay"))

#now rerun the PCA with just those groups
mcglPCA(subset.genind)
























### genind sample IDs. Ask Ben what the SRR numbers are. We're expecting the ##-###### format
name_df <-indNames(snps.genind.filt) %>% 
  as_tibble() %>% 
  separate(value, into = c("indiv", NA), sep = "_")

matching_indices <- match(indNames(snps.genind.filt), key$Run)


indNames(snps.genind.filt)[!is.na(matching_indices)] <- key$Sample.Name[matching_indices[!is.na(matching_indices)]]


### don't currently match the original IDs in
pop$indiv
