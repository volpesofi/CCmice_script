# an R script to perform QTL analysis with R/qtl2
# Authors: Anna Sofia Tascini

###### load libraries #######
library(plyr)
library("qtl2")
library("qtl2convert")
library("VariantAnnotation")
library("snpStats")
library("rgr")
library("broman")
library("readxl")
library("data.table")
library(patchwork)
library(cowplot)
library("optparse")
library(stringr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridisLite)
library("viridis")  
library(forcats)
library(ggplot2)
library(dplyr)

###### set wd #######
setwd("/beegfs/scratch/ric.cosr/ric.bragonzi/BragonziA_1641_WGS_mouse")


# --------- R qtl2 ------------
print("read cross")
FC_dataset = read_cross2("control.yaml", 
                         quiet = FALSE)

gmap <- insert_pseudomarkers(FC_dataset$gmap, step=0.1, stepwidth="max")
pmap <- interp_map(gmap, FC_dataset$gmap, FC_dataset $pmap)
ncore = 4

pr <- calc_genoprob(FC_dataset,
                    gmap,
                    map_function="c-f",
                    quiet = FALSE,
                    error_prob=0.002,
                    cores = ncore)

print("eval allele prob")
pr8 <- genoprob_to_alleleprob(probs = pr, 
                              quiet = FALSE, 
                              cores = ncore)

#### plotting ####
require(ggplot2)
require(reshape2)
print("Let's plot!")



# ------ plot -------
ylabel = "p"
correct = TRUE
  
plot_dir = "data/Mosaic_plot"
dir.create(plot_dir, recursive = T)
mouseID = row.names(FC_dataset$cross_info)[c(
  grep(pattern = "CC006_", row.names(FC_dataset$cross_info)),
  grep(pattern = "CC037_", row.names(FC_dataset$cross_info)))
]

chr_6_plots = list()
df_list = list()
correct_HETregion = TRUE

for (mouseID_i in 1:length(mouseID)) {
  mouseID_name = mouseID[mouseID_i]
  mouseID_name_mod = str_remove(string = mouseID_name, pattern = "\\/")
  print(paste("processing mouse", mouseID[mouseID_i]))
  plot = list()
  plot2 = list()
  XL = 200
  for (i in c(1:19, "X")) {
    chr_i = paste("chr", i, sep="")
    prob_m = pr8[[i]][mouseID_name,,]
    freqs = round(prob_m, digits = 5)
    freqdf <- as.data.frame(t(freqs))
    freqdf$pos = pmap[[i]]
    logodf <- data.frame(AJ=freqdf$A, B6=freqdf$B,
                         `129`=freqdf$C, NOD=freqdf$D,
                         NZO=freqdf$E, CAST=freqdf$F,
                         PWK=freqdf$G, WSB=freqdf$H,
                         pos=freqdf$pos, check.names = F)
    lmf <- melt(logodf, id.var='pos')
    colnames(lmf) <- c("pos", "founders", "value")
    lmf$r_value = round(lmf$value*2, digits = 0)/2
    
    
    plot[[chr_i]] = ggplot(data=lmf, aes(x=as.numeric(as.character(pos)), y=value))  +
      ylab(ylabel) + xlab(chr_i) +
      geom_bar(aes(fill=founders,
                   color = founders),
               position='stack',
               stat='identity', alpha=1) +
      scale_fill_manual(values = CCorigcolors) +
      scale_color_manual(values = CCorigcolors) +
      scale_x_continuous(limits = c(0,XL)) +
      geom_hline(yintercept=0.5, color = "white", size = 1) +
      theme_linedraw()
    
    plot2[[chr_i]] = ggplot(data=lmf, aes(x=as.numeric(as.character(pos)), y=r_value))  +
      ylab(ylabel) + xlab(chr_i) +
      ylim(c(0,1)) +
      geom_bar(aes(fill=founders,
                   color = founders),
               position='stack',
               stat='identity', alpha=1) +
      scale_fill_manual(values = CCorigcolors) +
      scale_color_manual(values = CCorigcolors) +
      scale_x_continuous(limits = c(0,XL)) +
      geom_hline(yintercept=0.5, color = "white", size = 1) +
      theme_linedraw() +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
    
    df_list[[mouseID_name_mod]][[chr_i]] = plot2[[chr_i]]$data
    
  }
  
  if (correct_HETregion) {
    
    if (mouseID_i == 1) {
      # correct het regions in CC006_S2_17Unc
      # from 0 to 32.64579 inherited from B6
      pos2MOD_1 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos <= 31.59669, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_1 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_1 & plot2[["chr6"]]$data$founders == "NZO","r_value"] = 0.5
      df_list[[mouseID_i]][["chr6"]] = plot2[["chr6"]]$data
    }

    
    if (mouseID_i == 2) {
      #correct het regions in CC037 based on het distribution
      
      # from 0 to 32.20757 het B6 + NOD
      
      pos2MOD_1 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos >=  0 & plot2[["chr6"]]$data$pos <= 32.2, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_1 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_1 & plot2[["chr6"]]$data$founders == "NOD","r_value"] = 0.5
      

      # 32.64579 to 34.81571 het B6 + WSB
      pos2MOD_2 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos > 32.2 & plot2[["chr6"]]$data$pos < 34.8157, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_2 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_2 & plot2[["chr6"]]$data$founders == "WSB","r_value"] = 0.5
      
      # 34.81571 to 37.75036 het B6 + CAST
      pos2MOD_3 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos >= 34.8157 & plot2[["chr6"]]$data$pos < 37.75036, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_3 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_3 & plot2[["chr6"]]$data$founders == "CAST","r_value"] = 0.5
      
      # 37.7503 to 39.94477 het B6 + NOD 
      pos2MOD_4 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos >= 37.7503 & plot2[["chr6"]]$data$pos <= 39.94477, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_4 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_4 & plot2[["chr6"]]$data$founders == "NOD","r_value"] = 0.5
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_4 & plot2[["chr6"]]$data$founders == "CAST","r_value"] = 0
      
      pos2MOD_5 = plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos >= 39.94477 & plot2[["chr6"]]$data$pos <= 49, "pos"]
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_5 & plot2[["chr6"]]$data$founders == "B6","r_value"] = 0
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_5 & plot2[["chr6"]]$data$founders == "NOD","r_value"] = 0
      plot2[["chr6"]]$data[plot2[["chr6"]]$data$pos %in% pos2MOD_5 & plot2[["chr6"]]$data$founders == "129","r_value"] = 1
      
      df_list[[mouseID_i]][["chr6"]] = plot2[["chr6"]]$data
    }
    

  }
  
  
  # save chr6 plot
  dir.create("data/Mosaic_plot/chr_6", recursive = T, showWarnings = F)
  c6p = chr_6_plots[[mouseID_name_mod]] = plot2[["chr6"]] + 
    xlim(c(0, 150)) + 
    geom_rect(mapping=aes(xmin=18, xmax=18.5, ymin=0, ymax=1), fill = "#CC9900",
              color = "#CC9900", alpha = 0.8) #+ legend_b
  ggsave(paste0("data/Mosaic_plot/chr_6/Haplotype_plot_", mouseID_name_mod, ".png"),
         plot = c6p,
         width = 20, height = 3)
  
  legend_b <- get_legend(
    plot[[chr_i]] +
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  
  patchwork <- plot[["chr1"]] / plot[["chr2"]] / plot[["chr3"]] / plot[["chr4"]] / plot[["chr5"]] / plot[["chr6"]] / plot[["chr7"]] /
    plot[["chr8"]] / plot[["chr9"]] / plot[["chr10"]] / plot[["chr11"]] / plot[["chr12"]] / plot[["chr13"]] /
    plot[["chr14"]] / plot[["chr15"]] / plot[["chr16"]] / plot[["chr17"]] / plot[["chr18"]] / plot[["chr19"]] / plot[["chrX"]]
  patchwork + plot_annotation(title = mouseID[mouseID_i],
                              theme = theme(plot.title = element_text(size = 18))) + plot_layout(guides = 'collect')
  
  ggsave(paste(plot_dir,"/mosaic_", mouseID_name_mod,".png", sep=''), height = 15, width = 10, bg = "white")
  
  patchwork <- plot2[["chr1"]] / plot2[["chr2"]] / plot2[["chr3"]] / plot2[["chr4"]] / plot2[["chr5"]] / plot2[["chr6"]] / plot2[["chr7"]] /
    plot2[["chr8"]] / plot2[["chr9"]] / plot2[["chr10"]] / plot2[["chr11"]] / plot2[["chr12"]] / plot2[["chr13"]] /
    plot2[["chr14"]] / plot2[["chr15"]] / plot2[["chr16"]] / plot2[["chr17"]] / plot2[["chr18"]] / plot2[["chr19"]] / plot2[["chrX"]]
  patchwork + plot_annotation(title = mouseID[mouseID_i],
                              theme = theme(plot.title = element_text(size = 18))) + plot_layout(guides = 'collect')
  
  ggsave(paste(plot_dir,"/mosaic_rounded_", mouseID_name_mod,".png", sep=''), height = 15, width = 10, bg = "white")
  
  v = data.frame()
  for (i in c(1:19, "X")) {
    chr_i = paste("chr", i, sep="")
    v = rbind(v,plot2[[chr_i]]$data[plot2[[chr_i]]$data$r_value>0,])
  }
  v$founders = as.character(v$founders)
  p = v %>%
    group_by(founders) %>%
    tally() %>%
    mutate(perc = n/sum(n)*100)  %>%
    ggplot(., aes(x = '', y = n, fill = founders)) +
    geom_bar(width = 1, stat = "identity", color = "black") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = CCorigcolors) +
    theme_void()
  print(p)
  tableF = as.data.frame(p$data)
  write.table(tableF,paste(plot_dir,"/percALL_table_", mouseID_name_mod,".txt", sep=''), quote = F, row.names = F, sep= "\t")
  ggsave(paste(plot_dir,"/percALL_", mouseID_name_mod,".png", sep=''), p, height = 15, width = 10, bg = "white")
  
}
