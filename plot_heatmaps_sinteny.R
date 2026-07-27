library(ggplot2)
library(dplyr)

species <- c("Species1","Species2","Species3","Species4","Species5")

all_list <- list()
outside_list <- list()
far_list <- list()

for (sp in species) {

  cat("Processing:", sp, "\n")

  sats <- read.table(paste0(sp, "_sats.bed"))
  colnames(sats) <- c("chr","start","end","family")
  sats$species <- sp

  outside <- read.table(paste0(sp, "_outside.bed"))
  colnames(outside) <- c("chr","start","end","family")
  outside$species <- sp

  dist_file <- paste0(sp, ".dist.txt")

  if (file.exists(dist_file)) {

    dist <- read.table(dist_file)

    if (ncol(dist) == 9) {
      colnames(dist) <- c("chr","start","end","family","chr2","start2","end2","block","dist")
    } else {
      colnames(dist) <- c("chr","start","end","family","chr2","start2","end2","dist")
    }

    far <- dist[dist$dist > 100000, ]
    far$species <- sp

  } else {
    far <- data.frame(chr=character(),start=integer(),end=integer(),family=character(),species=character())
  }

  all_list[[sp]] <- sats
  outside_list[[sp]] <- outside
  far_list[[sp]] <- far
}

all_df <- bind_rows(all_list)
outside_df <- bind_rows(outside_list)
far_df <- bind_rows(far_list)

total <- all_df %>% group_by(species, family) %>% summarise(total=n(), .groups="drop")
outside <- outside_df %>% group_by(species, family) %>% summarise(outside=n(), .groups="drop")
far <- far_df %>% group_by(species, family) %>% summarise(far=n(), .groups="drop")

summary_df <- merge(total, outside, all=TRUE)
summary_df <- merge(summary_df, far, all=TRUE)
summary_df[is.na(summary_df)] <- 0
summary_df$perc_outside <- (summary_df$outside / summary_df$total) * 100
summary_df$perc_far <- (summary_df$far / summary_df$total) * 100

summary_df$species <- factor(summary_df$species, levels=species)

pdf("heatmap_outside.pdf",6,5)
print(
ggplot(summary_df,aes(x=species,y=family,fill=perc_outside))+
geom_tile()+scale_fill_gradient(low="white",high="red")+theme_minimal()+
labs(title="% satDNA outside synteny"))
dev.off()

pdf("heatmap_far.pdf",6,5)
print(
ggplot(summary_df,aes(x=species,y=family,fill=perc_far))+
geom_tile()+scale_fill_gradient(low="white",high="blue")+theme_minimal()+
labs(title="% satDNA far (>100kb)"))
dev.off()

cat("Done! Plots generated.\n")
