### Challenge 3 version by Fran. Jul 2021

LRGASP_calculations_challenge3 <- function (NAME, class.file, junc.file, out.dir, platform, functions.dir, bam) {
  # Get functions and spike-ins IDs
  setwd(functions.dir)
  source("LRGASP_functions.R")
  sirv_list=read.table("SIRVs_ids.txt", header = F)$V1

  # identify files in directory
  cat("Evaluation script has being run.\nData used for ", NAME, " pipeline are \n", class.file , "\n", junc.file , "\n")
  sqanti_data=read.table(class.file , sep = "\t", as.is = T, header = T)
  sqanti_data.junc=read.table(junc.file, sep = "\t", as.is = T, header = T)
  if(all(is.na(sqanti_data$iso_exp))){
    sqanti_data$iso_exp <- 0
  }
  # change names of structural categories
  cat_levels <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
  cat_labels <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  sqanti_data$structural_category = factor(sqanti_data$structural_category,
                                          labels = cat_labels,
                                          levels = cat_levels,
                                          ordered=TRUE)

  ### Add LRGASP_id
  
  iso_tags <- isoformTags(sqanti_data.junc)
  sqanti_data <- merge(sqanti_data, iso_tags, by="isoform" , all.x=TRUE)
  sqanti_data$LRGASP_id <- apply(sqanti_data, 1, monoexon_tag)
  sqanti_data <- addSC(sqanti_data)

 ### rewrite initial table with LRGASP_id
  
  write.table(sqanti_data, class.file, quote=F, sep = "\t", row.names = FALSE)

# separate spike-ins and isoforms
  sirv_data=sqanti_data[grep("SIRV",sqanti_data$chrom),]
  sirv_data.junc=sqanti_data.junc[grep("SIRV",sqanti_data.junc$chrom),]

### remove SIRV and ERCC transcripts from sqanti data
  sqanti_data=sqanti_data[grep("SIRV",sqanti_data$chrom, invert=T),]
  sqanti_data.junc=sqanti_data.junc[grep("SIRV",sqanti_data.junc$chrom, invert=T),]
  
  ### Evaluation of isoforms
  #####################
  print ("Evaluation without a reference transcriptome")
  bamFile <- BamFile(bam)
  align_file <- scanBam(bamFile)
  flags_t <- as.data.frame(table(align_file[[1]]$flag))
  num_reads_map_SIRV <- length(grep(align_file[[1]]$rname, pattern = "SIRV"))
  rownames(flags_t) <- flags_t$Var1
  # mapping rate is calculated without taking into account reads mapping to SIRV sequences
  num_isoforms = sum(flags_t$Freq)-num_reads_map_SIRV
  mapping_isoforms=sum(flags_t$Freq)-num_reads_map_SIRV-flags_t["4","Freq"]
  mapping_rate=(1-(flags_t["4","Freq"]/num_isoforms))*100
  
  # length isoforms
  average_length=mean(sqanti_data$length)
  sd_length=sd(sqanti_data$length)
  
  # coding potential
  coding_isoforms=length(which(sqanti_data$coding=="coding"))
  perc_coding=(coding_isoforms/num_isoforms)*100
  
  # Full Illumina SJ Support
  full_Illumina_SJ_support=length(which(sqanti_data$min_cov>0))
  perc_full_Illumina_SJ_support=(full_Illumina_SJ_support/num_isoforms)*100
  
  # Non canonical Isoforms
  non_canonical_isoforms=length(which(sqanti_data$all_canonical=="non_canonical"))
  perc_non_canonical_isoforms=(non_canonical_isoforms/num_isoforms)*100
  
  # Intra-priming
  intrapriming_isoforms=length(which(sqanti_data$perc_A_downstream_TTS>=60))
  perc_intraprimimg_isoforms=(intrapriming_isoforms/num_isoforms)*100
  
  # RT-switching
  RTS_isoforms=length(which(sqanti_data$RTS_stage=="TRUE"))
  perc_RTS_isoforms=(RTS_isoforms/num_isoforms)*100
  
  # Illumina coverage
  select_cols=c("chrom", "strand", "genomic_start_coord", "genomic_end_coord", "junction_category", "start_site_category", 
                "end_site_category", "splice_site", "canonical", "sample_with_cov", "total_coverage")
  unique_SJ=unique(sqanti_data.junc[,select_cols])
  SJ_with_cov=unique_SJ[which(unique_SJ$total_coverage>0),]
  num_SJ_with_cov=length(SJ_with_cov$strand)
  perc_SJ_with_cov=(num_SJ_with_cov/length(unique_SJ$strand))*100
  
  # SJ canonical
  SJ_non_canonical=unique_SJ[which(unique_SJ$canonical=="non_canonical"),]
  num_SJ_non_canonical=length(SJ_non_canonical$strand)
  perc_SJ_non_canonical=(num_SJ_non_canonical/length(unique_SJ$strand))*100
  
  
  
  
  # Write out results
  a.non_model_results=data.frame(row.names = c("Number of transcripts", "Mapping transcripts", "Average length", "Transcripts with coding potential",
                                               "Transcripts with Full Illumina SJ Support", "Non-canonical transcripts",
                                               "Transcripts with possible intra-priming", "Transcripts with possible RT-switching",
                                               "Splice Junctions with short-read coverage", "Non-canonical Splice Junctions"))
  
  a.non_model_results[,"Absolute value"]="-"
  a.non_model_results[,"Relative value (%)"]="-"
  
  a.non_model_results["Number of transcripts","Absolute value"]=as.integer(num_isoforms)
  a.non_model_results["Mapping transcripts","Absolute value"]=as.integer(mapping_isoforms)
  a.non_model_results["Mapping transcripts","Relative value (%)"]=round(mapping_rate, digits=2)
  a.non_model_results["Average length","Absolute value"]=round(average_length, digits=2)
  a.non_model_results["Transcripts with coding potential","Absolute value"]=as.integer(coding_isoforms)
  a.non_model_results["Transcripts with coding potential", "Relative value (%)"]=round(perc_coding, digits=2)
  a.non_model_results["Transcripts with Full Illumina SJ Support","Absolute value"]=as.integer(full_Illumina_SJ_support)
  a.non_model_results["Transcripts with Full Illumina SJ Support", "Relative value (%)"]=round(perc_full_Illumina_SJ_support, digits=2)
  a.non_model_results["Non-canonical transcripts","Absolute value"]=as.integer(non_canonical_isoforms)
  a.non_model_results["Non-canonical transcripts", "Relative value (%)"]=round(perc_non_canonical_isoforms, digits=2)
  a.non_model_results["Transcripts with possible intra-priming","Absolute value"]=as.integer(intrapriming_isoforms)
  a.non_model_results["Transcripts with possible intra-priming", "Relative value (%)"]=round(perc_intraprimimg_isoforms, digits=2)
  a.non_model_results["Transcripts with possible RT-switching","Absolute value"]=as.integer(RTS_isoforms)
  a.non_model_results["Transcripts with possible RT-switching", "Relative value (%)"]=round(perc_RTS_isoforms, digits=2)
  
  a.non_model_results["Splice Junctions with short-read coverage", "Absolute value"]=as.integer(num_SJ_with_cov)
  a.non_model_results["Splice Junctions with short-read coverage", "Relative value (%)"]=round(perc_SJ_with_cov, digits=2)
  a.non_model_results["Non-canonical Splice Junctions", "Absolute value"]=as.integer(num_SJ_non_canonical)
  a.non_model_results["Non-canonical Splice Junctions", "Relative value (%)"]=round(perc_SJ_non_canonical, digits=2)
  
  ### Evaluation of SIRVs
  ##############################################
  print ("SIRVs evaluation")
  sirv_data$TP=apply(sirv_data,1,TP_function)
  SIRVs_transcripts=as.integer(length(sirv_data$isoform))
  SIRVs_called=intersect(sirv_data[which(sirv_data$structural_category=="FSM" & 
                                           sirv_data$TP==TRUE),"associated_transcript"],
                         sirv_list)
  TP=length(SIRVs_called)
  RM_isoforms=sirv_data[which(sirv_data$structural_category=="FSM" & 
                 sirv_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  SIRVs_transcripts_incomplete=sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM") & 
                                                 sirv_data$TP==FALSE),"isoform"]
  
  SIRVs_called_wrong_ends=setdiff(intersect(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"],
                                    sirv_list),SIRVs_called)
  PTP=length(SIRVs_called_wrong_ends)
  SIRVs_not_detected=setdiff(sirv_list,sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"])
  FN=length(SIRVs_not_detected)
  FP_sirvs_detected=sirv_data[-which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"]
  FP=length(FP_sirvs_detected)
  
  SIRVs_redundancy=length(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"])/length(unique(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"]))
  

  # Write out results
  b.SIRVs_results=data.frame(row.names = c("SIRV transcripts", "True Positive detections (TP)", "SIRV transcripts associated to TP (Reference Match)",
                                           "Partial True Positive detections (PTP)", "SIRV transcripts associated to PTP",
                                           "False Negative (FN)", "False Positive (FP)", 
                                           "Sensitivity", "Precision",
                                           "Non Redundant Precision","Positive Detection Rate",
                                           "False Discovery Rate","Redundancy"))
  b.SIRVs_results[,"Value"]="-"
  b.SIRVs_results["SIRV transcripts","Value"]=SIRVs_transcripts
  b.SIRVs_results["True Positive detections (TP)","Value"]=as.integer(TP)
  b.SIRVs_results["SIRV transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  b.SIRVs_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  b.SIRVs_results["SIRV transcripts associated to PTP","Value"]=as.integer(length(SIRVs_transcripts_incomplete))
  b.SIRVs_results["False Negative (FN)","Value"]=as.integer(FN)
  b.SIRVs_results["False Positive (FP)","Value"]=as.integer(FP)
  b.SIRVs_results["Sensitivity","Value"]=round(TP/length(sirv_list), digits = 2)
  b.SIRVs_results["Precision","Value"]=round(RM/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Non Redundant Precision","Value"]=round(TP/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Positive Detection Rate", "Value"]=round(length(unique(c(SIRVs_called,SIRVs_called_wrong_ends)))/length(sirv_list), digits = 2)
  b.SIRVs_results["False Discovery Rate","Value"]=round((SIRVs_transcripts - RM)/SIRVs_transcripts, digits = 2)
  b.SIRVs_results["Redundancy","Value"]=round(SIRVs_redundancy, digits = 2)
                              
                  
  
  ####Create a list with all results and save all
  ###############################################
  
  files <- ls(pattern = "_results")
  all.results <- list()
  for ( i in 1: length(files) ) {
    all.results[[i]] <- eval(parse(text = files[i]))
  }
  setwd(out.dir)
  names(all.results) <- c("Transcriptome without reference", "SIRV") 
  
  save(all.results , file = paste(NAME, "_results.RData", sep = ''))
  save(sqanti_data, file=paste(NAME, "_classification.RData", sep = ''))
  save(sqanti_data.junc, file=paste(NAME, "_junctions.RData", sep = ''))
  save(sirv_data, file=paste(NAME, "_SIRVs_class.RData", sep=''))
  save(sirv_data.junc, file=paste(NAME, "_SIRVs_junc.RData", sep=''))

}

