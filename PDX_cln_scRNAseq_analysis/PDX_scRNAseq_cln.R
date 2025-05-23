###
# PDX scRNA-seq: Analysis of Long Nguyen's scRNA-seq barcoded data
# 
###

library(plyr, quietly = T)
library(RColorBrewer, quietly = T)
library(mclust, quietly = T)
library(doMC, quietly = T)
library(tidyverse, quietly = T)
library(tgstat, quietly = T)
library(Matrix, quietly = T)
library(pheatmap, quietly = T)
library(dplyr, quietly = T)
library(metacell, quietly = T)

source("LN_mc_utils.R")

# change working directory if needed to seperate from code directory 
wd = "."

rl(paste0(wd, "/scdb"), paste0(wd, "/figs"), config_fn="config/pdx_Long.yaml", force_init=F)

# Colour defs ----
md_cols = list(Model = c(STG139="red", STG201="orange", NKI250="magenta", AB040="blue", IC07="darkgreen"),
               Subtype = setNames(c('blue', 'red'), c('ER-positive', 'TN')),
               Replicate = c(P1="#4292C6", P2="#08306B", S1="#FEE391", S2="#EC7014", S3="#662506"),
               Passage = c(X4="#41AB5D", X5="#00441B"),
               CloneType = c(Propagating='red', Transient='darkgray', Emerging='blue'))

big_cln_cols = c(STG139_cl1 = "#984EA3", 
                STG201_cl6 = "#377EB8", 
                AB040_cl40=  "#E41A1C", 
                AB040_cl76= "#FF7F00" , 
                AB040_cl258= "#FFFF33", 
                IC07_cl273 = "#4DAF4A")

model2subtype = c(STG139="TN", STG201="TN", NKI250="TN", IC07="ER-positive", AB040="ER-positive")

LN_get_sample_info = function(report_ifn="report_v2.csv", selected_columns=c('SampleID', 'Sample.name', 'Run_name', 'Loaded_viable_cells',	'Model', 'Replicate',	'CellID_suffix', 'Mouse_ID', 'Passage'))
{
  read.csv(report_ifn) %>%
    dplyr::select(all_of(selected_columns))
  
}

#
# Aggregate samples count matrices into a single mat ----
# Change runs_dir to base dir of 10x processed runs
LN_build_mat_v2 = function(mat_id, runs_dir = "/mnt/scratcha/cclab/scRNAseq/runs", samples_info_ifn="report_v2.csv")
{
  prev_scdb = curr_scdb = .scdb_base
  
  sinfo = LN_get_sample_info(report_ifn = samples_info_ifn)
  
  aggr_csv = NULL
  
  mat = NULL
  md = NULL
  
  for (i in 1:nrow(sinfo)) {
    nm = sinfo[i, 'Sample.name']
    idir = sprintf("%s/%s/scdb", runs_dir, sinfo[i, 'Run_name'])
    if (idir != curr_scdb) {
      scdb_init(idir, T)
      curr_scdb = .scdb_base
    }
    
    cmat = scdb_mat(nm)
    if (!is.null(cmat) && cmat@ncells > 0) {
      
      cells = gsub("-.*", paste0("-", sinfo[i, 'CellID_suffix']), cmat@cells)
      
      suff = unique(gsub(".*-", "", cmat@cells))
      new_suff = unique(gsub(".*-", "", cells))
      message(glue("Adding mat {i} {nm} with suff {suff} -> {new_suff}..."))
      
      colnames(cmat@mat) = cells
      rownames(cmat@cell_metadata) = cells
      
      cmd = cbind(cmat@cell_metadata, sinfo[i,])
      
      cmat@mat = cmat@mat[grepl("GRCh38_", rownames(cmat@mat)), ]
      rownames(cmat@mat) = gsub("GRCh38_", "", rownames(cmat@mat))
      
      if (!is.null(mat)) {
        stopifnot(all(rownames(mat) == rownames(cmat@mat)))
      }
      
      mat = cbind(mat, cmat@mat)
      md = rbind(md, cmd)
      
    } else {
      message(sprintf("%2d: skipping %s", i, nm))
    }
  }
  
  md$batch_set_id = md$cDNA_prep_date
  md$seq_batch_id = md$Run_name
  
  scdb_init(prev_scdb, T)
  message(sprintf("==== Creating %s mat with %d cells ====", mat_id, ncol(mat)))
  
  scdb_add_mat(paste0(mat_id, "_full"),  tgScMat(mat, "umi", md))
  
  # filter by umis + blacklisted cells
  bl_cells = read.csv("blacklist_cells.csv", header=F)$V1
  valid_cells = setdiff(rownames(md)[md$valid_by_umis], bl_cells)
  
  mat = mat[, valid_cells]
  md = md[valid_cells, ]
  scdb_add_mat(mat_id,  tgScMat(mat, "umi", md))
  
}

# Generare and annotate metacell model -----
LN_pipe_v2 = function(rebuild=F)
{
  set.seed(42)
  
  ds_nm = "pdx_LN_v2" 
  filt_id = sprintf("%s_filt", ds_nm)
  lateral_gset_id = sprintf("%s_lateral", ds_nm)
  
  # Use all cells (already filtered cells by QC)
  common_pipe(ds_nm=ds_nm, LN_build_mat_v2, max_f_mit=1, min_umis_pre_gene_ignore=0, min_umis_post_gene_ignore=0, gsets_to_filter=c("mito", "ncrna"),
              gset_cor_thresh=0.1,  
              T_lfc = 3000, cgraph_downsamp=T, 
              bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
              mc_K=30, min_mc_size=20, mc_alpha=2, sample_field="SampleID", 
              color_by_conf=F, rebuild=rebuild, out_base=NULL, specie='human',split_mc_with_dbscan_and_filt_outliers=F)
  
  # Annotate mcs by normal epithelial modules from Pal 2021
  sinfo = LN_get_sample_info()
  models = sort(unique(sinfo$Model))
  model_ann = data.frame(row.names=models, Model=models)
  
  pdx_annotate_mcs_by_normal_signatures(mc_id = filt_id, samp_ann = model_ann, samp_ann_cols = md_cols, analysis_name = ds_nm,
                                        pal2021_epi_marks_method="tot_group_umis", 
                                        rebuild=rebuild)
}


##
#
pdx_annotate_mcs_by_normal_signatures = function(mc_id, samp_ann, samp_ann_cols, analysis_name,
                                                 sample_md_col = "SampleID", annotate_mc_by_col='Model',
                                                 mat_id = mc_id, graph_id = mat_id, raw_mat_id = gsub("filt", "raw", mat_id),
                                                 mc_col_by_cl_id = sprintf("%s_colByCluster_%s", mc_id, analysis_name),
                                                 cl_cols = setNames(c('#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#444444', '#CCCCCC'), c('Basal', 'LP', 'ML', 'Basal-LP', 'Basal-ML', 'LP-ML', 'Basal-LP-ML', 'None')),
                                                 rebuild = F)
{
  curr_fdir = .scfigs_base
  odir = scfigs_dir(mc_id, sprintf("%s_deBy%s_qForPosSig%.2f", analysis_name, pal2021_epi_marks_method, sig_min_q))
  scfigs_init(odir)
  
  # load data
  mat = scdb_mat(mat_id)
  mc = scdb_mc(mc_id)
  
  lfp = log2(mc@mc_fp)
  col2grp = get_mc_col2group(mc)
  
  # prepaing metadata
  mixed_ann = t(as.data.frame(setNames(rep(NA, ncol(samp_ann)), colnames(samp_ann))))
  rownames(mixed_ann) = 'Mixed'
  mc_type_ann = rbind(samp_ann, mixed_ann)

  mc_samp = table(mc@mc, mat@cell_metadata[names(mc@mc), annotate_mc_by_col])
  mc_samp_n = mc_samp / rowSums(mc_samp)
  mc2samp = ifelse(apply(mc_samp_n, 1, max) >= 0.9, colnames(mc_samp_n)[apply(mc_samp_n, 1, which.max)], 'Mixed')
  mc_type_df = mc_type_ann[mc2samp, , drop=F]
  rownames(mc_type_df) = colnames(lfp)
  
  # Annotate mcs by Pal2021 epihelial markers (and lateral) ---
  pal_sig_gset = scdb_gset(mc_col_by_cl_id)@gene_set
  supmc2marks = lapply(split(names(pal_sig_gset), pal_sig_gset), intersect, y=rownames(lfp))
  
  # annotate by mean lfp (mc level)
  mc_avg_epi_sig = t(sapply(supmc2marks, function(v) colMeans(lfp[intersect(rownames(lfp), v), ])))
  mc_max_sig = rownames(mc_avg_epi_sig)[apply(scale(t(mc_avg_epi_sig)), 1, which.max)]
  mc_sig_q = sapply(rownames(mc_avg_epi_sig), function(nm) as.numeric(quantile(mc_avg_epi_sig[nm, mc_max_sig == nm], sig_min_q)))
  mc_cls = apply(sapply(sort(rownames(mc_avg_epi_sig)), function(sig) { mc_avg_epi_sig[sig, ] >= mc_sig_q[sig]} ), 1, function(v) paste0(sort(rownames(mc_avg_epi_sig))[v], collapse="-"))
  mc_type_df$cl = ifelse(mc_cls == "", "None", mc_cls)
  
  mc_ord = NULL
  cl_sizes = NULL
  for (cl in intersect(names(cl_cols), unique(mc_type_df$cl))) { 
    ind = mc_type_df$cl == cl
    hc = hclust(dist(t(mc_avg_epi_sig[,ind])), method='ward.D2')
    mc_ord = c(mc_ord, colnames(mc_avg_epi_sig)[ind][hc$order])
    cl_sizes = c(cl_sizes, sum(ind))
  }
  
  ann_cols = samp_ann_cols
  ann_cols[['cl']] = cl_cols
  zlim = 4; ch = 20; cw=1;

  png(scfigs_fn(mc_id, "avg_normal_Epithelial_by_mc_cutoffs"), ncol(mc_avg_epi_sig) * cw + 300, nrow(mc_avg_epi_sig) * ch + 400)
  pheatmap(pmin(pmax(t(scale(t(mc_avg_epi_sig[c('Basal', 'LP', 'ML'), mc_ord]))), -zlim), zlim), breaks=seq(-zlim, zlim, len=101), color=colorRampPalette(c("#313695", "white", "#A50026"))(100), cluster_cols=F, gaps_col=cumsum(cl_sizes), annotation_col=mc_type_df, annotation_colors=ann_cols, cluster_rows=F, treeheight_col=20, cellwidth=cw, cellheight=ch, show_colnames=F)
  dev.off()
  
  # signature by cluster (is it only a matter of cutoffs? e.g. basal-lm is medium in both basal and lm compared to pure basal and pure lm?)
  sig_avg_df = reshape2::melt(mc_avg_epi_sig, varnames=c('sig', 'mc'), value.name = 'mean') %>%
    mutate(sig = factor(sig, levels=c('Basal', 'LP', 'ML')))
  
  sig_avg_df$cl = mc_type_df[sig_avg_df$mc, 'cl']
  
  p = ggplot(sig_avg_df, aes(x=cl, y=mean, fill=cl)) + geom_boxplot(notch=T, outlier.shape=19, outlier.size=0.2, show.legend=F) + facet_wrap(~sig, nrow=1, scales="free_y") + labs(x="", y="Mean signature") + scale_fill_manual(values=cl_cols) + scale_x_discrete(guide=guide_axis(angle=30))
  ggsave(scfigs_fn(mc_id, "mean_epithelial_signatures_by_cluster"), p, width=10, height=4)
  
  # annotate mcs by clusters
  mc_col_by_cl = mc
  mc_col_by_cl@colors = cl_cols[mc_type_df$cl]
  mc_col_by_cl@color_key = data.frame(row.names = seq_along(cl_cols), gene="", group=names(cl_cols), color=cl_cols)
  if (!scdb_obj_exists("mc", mc_col_by_cl_id) || rebuild) {
    scdb_add_mc(mc_col_by_cl_id, mc_col_by_cl)
  }
  mcell_mc2d_force_knn(mc_col_by_cl_id, mc_col_by_cl_id, graph_id)
  mcell_mc2d_plot(mc_col_by_cl_id, plot_box_and_axes=F, show_mc_ids=F)
  
  col2type = get_mc_col2group(mc_col_by_cl)

  if (!scdb_obj_exists("mat", mc_col_by_cl_id) || rebuild) {
    mat_col_by_cl = mat
    mat_col_by_cl@cell_metadata$supmc = NA
    mat_col_by_cl@cell_metadata[names(mc_col_by_cl@mc), 'supmc'] = col2type[mc_col_by_cl@colors[mc_col_by_cl@mc]]
    
    scdb_add_mat(mc_col_by_cl_id, mat_col_by_cl)
  }
  
  scfigs_init(curr_fdir)
  
  invisible(col2type[mc_col_by_cl@colors])
}

# Analysis ----

#
# Simulate clones competetion and recovery by actual number of initiating cells and sequenced reads 
##
LN_CIC_analysis = function(info_ifn = "CIC_info.csv", n_cores = 8, model_plot_ord = c('STG139', 'STG139M', 'STG201', 'AB040', 'IC07'),
                           f_init_death=0, f_death=0, f_double=1, n_sim=50, rebuild=F)
{
  x = read.csv(info_ifn) %>%
    rename(Sample_name = 'Model') %>%
    mutate(Model = gsub("_.*", "", Sample_name))

  x_f = group_by(x, Model) %>%
    summarise(n=n()) %>%
    filter(n > 4) %>%
    inner_join(x)
  
  ofn_pref = glue("{.scfigs_base}/CIC_fInitDeath_{f_init_death}_fDeath_{f_death}_fDouble_{f_double}_{n_sim}_iters")
  
  # simulation function
  inner_sim_expansion = function(j, f_init_death, f_death, f_double, target_cells, n_sim) {
    n_imp = x_f$Cells.implanted[j]
    f_gfp = x_f$X.GFP[j]/100
    n_reads = x_f$Total.reads[j]
    
    n_gfp_clones = c()
    
    message(glue("{j}: {x_f$Sample_name[j]} ({f_init_death}, {f_death}, {f_double}, {n_imp}, {f_gfp}, {n_reads})"))
    
    for (i in 1:n_sim) {
      v = 1:n_imp
      stopifnot(n_imp < target_cells)
      
      # Initial death (simulating failure to implant)
      v = sample(v, size = floor(length(v) * (1 - f_init_death)), replace = F)
      
      # Death+expand cycles until reaching target number of cells
      while (length(v) <= target_cells) {
        #message(glue("Iter {i}: {length(v)} cells..."))
        # lose cells that die
        v = sample(v, size = floor(length(v) * (1 - f_death)), replace = F)
        
        # cells dividing (doubling selected f_double cells)
        v = c(v, sample(v, size = floor(length(v) * f_double), replace = F))
      }
      # simulate selectiong of reads (with replacement, mimicking PCR amplification)
      seq_v = sample(v, n_reads, replace=T)
      
      # Filtering for cells with gfp
      seq_v = seq_v[seq_v <= round(n_imp * f_gfp)]
      
      # count distinct clones
      n_gfp_clones = c(n_gfp_clones, length(unique(seq_v)))
      #message(glue("Iter {i}: {length(unique(seq_v))} GFP clones"))
    }
    return(data.frame(mean_n_clones=mean(n_gfp_clones), sd_n_clones=sd(n_gfp_clones)))
  }
  
  # run the simulations
  res_fn = paste0(ofn_pref, "_res.csv")
  if (rebuild || !file.exists(res_fn)) {
    doMC::registerDoMC(cores = n_cores)
    res = alply(1:nrow(x_f), 1, inner_sim_expansion, f_init_death = f_init_death, f_death = f_death, f_double = f_double, target_cells = mean(x_f$Estimated.total.cells.per.tumour), n_sim = n_sim, .parallel = n_cores > 1)
    stopifnot(all(sapply(res, typeof) == 'list'))
    
    x_f_r = cbind(x_f, do.call('rbind', res)) %>%
    mutate(sim_CIC_freq = mean_n_clones / (Cells.implanted * (X.GFP / 100)))
  
    write.csv(x_f_r, res_fn, quote=F, row.names = F)
  }
  x_f_r = read.csv(res_fn) %>%
    mutate(n_clones_z = (Clones.detected - mean_n_clones) / sqrt(mean_n_clones),
           model_ord = factor(Model, levels=model_plot_ord))
  
  gg_df = dplyr::select(x_f_r, Cells.implanted, model_ord, sim_CIC_freq, CIC.frequency) %>% 
    rename(Observed = 'CIC.frequency', Simulated = 'sim_CIC_freq') %>%
    pivot_longer(cols=3:4, names_to='type', values_to='CIC_freq')
    
  
  p = ggplot(gg_df, aes(x=log10(Cells.implanted), y = log10(CIC_freq), color=type)) + geom_point(shape=19) + geom_smooth(method = 'lm', se = T) + facet_wrap(~model_ord, nrow=1) + labs(title = glue("Init %death: {f_init_death}, %death: {f_death}, %double: {f_double}"))
  ggsave(paste0(ofn_pref, "_simulated_cic_vs_nImp.png"), p, width=9, height=3)

  slopes = matrix(0, nrow=length(unique(gg_df$type)), ncol=length(unique(gg_df$model_ord)), dimnames=list(unique(gg_df$type), levels(gg_df$model_ord)))
  for (c_model in unique(gg_df$model_ord)) {
    for (c_type in unique(gg_df$type)) {
      gg_df_f = filter(gg_df, model_ord == c_model & type == c_type)
      slopes[c_type, c_model] = lm(log10(CIC_freq) ~ log10(Cells.implanted), gg_df_f)$coef[2]
    }
  }
  write.csv(slopes, paste0(ofn_pref, "_slopes.csv"), quote=F)
  
  p = ggplot(x_f_r, aes(x=log10(Cells.implanted), y = n_clones_z)) + geom_point(shape=19) + geom_smooth(method = 'lm', se = T) +  facet_wrap(~model_ord, nrow=1) + labs(title = glue("Init %death: {f_init_death}, %death: {f_death}, %double: {f_double}"), ylab="#clones (z-score)")
  ggsave(paste0(ofn_pref, "_simulated_cic_vs_nImp_zScore.png"), p, width=9, height=3)
 
}

# 
# Per large clone, create a metacell model (mc + mc2d) and order metacells by distance from metacells with primary cells and find genes(/modules) correlated with that ordering
###
LN_clones_dist_from_primary_analysis = function(mat_id = "pdx_LN_v2_filt_colByCluster_pdx_LN_v2",
                                                bc_info_csv = "LN_bc_info.csv",
                                                foc_bcs = c(STG139_cl1 = 'GFP_ACTGATATCGAGATCGAAACTGGTAAAACGT', 
                                                            STG201_cl6 = 'GFP_ACTGTAATCTAGATCGAAATAGGTAAAACCC', 
                                                            AB040_cl40='GFP_ACTGTAATCAGGATCCAAACTGGTGTAACGT', 
                                                            AB040_cl76='GFP_ACTGGGATCACGATGCAAAAGGGTGTAACAC', 
                                                            AB040_cl258='GFP_ACTGAAATCAAGATGGAAAGAGGTCCAACAT', 
                                                            IC07_cl273='GFP_ACTGCCATCCTGATCCAAACCGGTTAAACAA'),
                                                lateral_gset_id = "pdx_LN_v2_lateral",
                                                min_cells_for_primary_mc = c(STG139_cl1 = 3, 
                                                            STG201_cl6 = 3, 
                                                            AB040_cl40=2, 
                                                            AB040_cl76=2, 
                                                            AB040_cl258=2, 
                                                            IC07_cl273=5),
                                                min_lfp_max = 0.7, disp_n_d_cor_gs=5, min_gm_mu_diff = 0.3, n_hc_sd_for_cutree = 3, 
                                                gsea_max_padj = 0.05, 
                                                tfs_ifn = "Lambert2018_TF_names_v_1.01.txt",
                                                epi_gset = "pdx_LN_v2_filt_colByCluster_pdx_LN_v2",  
                                                supmc_cols = setNames(c('#FF0000', '#00FF00', '#0000FF', '#FFFF00', '#FF00FF', '#00FFFF', '#444444', '#CCCCCC'), c('Basal', 'LP', 'ML', 'Basal-LP', 'Basal-ML', 'LP-ML', 'Basal-LP-ML', 'None')),
                                                cyclone_cols = c(G1="#8DD3C7", S="#FFFFB3", G2M="#BEBADA"),
                                                rebuild=F
                                                )
{
  sinfo = LN_get_sample_info()
  
  mat = scdb_mat(mat_id)
  bc_info = fread(bc_info_csv, sep=",")
  
  md = mat@cell_metadata %>% 
    rownames_to_column(var='cellid') %>%
    left_join(bc_info) %>%
    filter(!is.na(var_seq))
  
  bc_info = bc_info %>% column_to_rownames(var='cellid')
  
  # generate a metacell model per clone
  for (foc_bc in names(foc_bcs)[1]) {
    message(glue("Processing clone {foc_bc}..."))
    foc_cells = rownames(bc_info)[bc_info$var_seq == foc_bcs[foc_bc]]
    
    c_id = paste(mat_id, foc_bc, sep="_")
    c_m = mat@mat[, foc_cells]
    c_md = cbind(mat@cell_metadata[foc_cells, ], bc_info[foc_cells, ])
    scdb_add_mat(c_id, tgScMat(c_m, "umi", c_md))
    
    if (rebuild || !scdb_obj_exists("mc", c_id)) {
      meta_mat2mc(c_id, lateral_gset_id=lateral_gset_id, 
                  split_mc_with_dbscan_and_filt_outliers=F, rebuild=rebuild) 
    }
    
    lat_genes = names(scdb_gset(lateral_gset_id)@gene_set)
    tgconfig::set_param("mcell_mc2d_max_confu_deg", 3, "metacell")
    
    if (rebuild || !scdb_obj_exists("mc2d", c_id)) {
      mcell_mc2d_force_knn(c_id, c_id, c_id,
                           ignore_mismatch=F, symmetrize=F, 
                           feats_gset = c_id, 
                           feats_exclude = lat_genes,
                           graph_parametric = T)
    }
  }
  
  tgconfig::set_param("mcell_mgraph_max_confu_deg", 3, "metacell")
  
  tfs = fread(tfs_ifn, header=F)$V1
  glob_gm_df = NULL
  glob_mu_gm_dist = NULL
  
  for (foc_bc in names(foc_bcs)[1]) {
    message(glue("Processing clone {foc_bc}..."))
    c_id = paste(mat_id, foc_bc, sep="_")
    
    mat = scdb_mat(c_id)
    mc = scdb_mc(c_id)
    mc2d = scdb_mc2d(c_id)
    mgraph = mc2d@graph
    
    n_primary_per_mc = tapply(grepl('^P', mat@cell_metadata[names(mc@mc), 'Replicate']), mc@mc, sum)
    mcs_with_primary = which(n_primary_per_mc >= min_cells_for_primary_mc[foc_bc])
    message(glue("{c_id}: {length(mcs_with_primary)} mcs with at least {min_cells_for_primary_mc[foc_bc]} cells"))
    
    order_mgraph_mcs_by_mc_dist = function(mgraph, src_mc, dists, symmetrise=T) {
      if (symmetrise) {
        mgraph1 = dplyr::select(mgraph, mc1, mc2) 
        mgraph2 = data.frame(mc1=mgraph$mc2, mc2=mgraph$mc1)
        mgraph = rbind(mgraph1, mgraph2) %>% distinct
      }
      
      dists[src_mc] = 0
      src_mcs = src_mc
      
      while(length(src_mcs) > 0 && any(is.na(dists))) {
        c_src_mc = src_mcs[1]
        #message(glue("c_src_mc: {c_src_mc} ({dists[c_src_mc]})\tTo fill: {sum(is.na(dists))}"))
        adj = filter(mgraph, mc1 == c_src_mc & (is.na(dists[mc2]) | dists[mc2] > dists[c_src_mc]))
      
        src_mcs = setdiff(src_mcs, c_src_mc)
        if (nrow(adj) > 0) {
          dists[adj$mc2] = ifelse(is.na(dists[adj$mc2]), dists[c_src_mc] + 1, pmin(dists[adj$mc2], dists[c_src_mc] + 1))
          #message(glue("Done with {c_src_mc}, adding: {paste0(adj$mc2, collapse=', ')}"))
          src_mcs = c(src_mcs, adj$mc2)
        }
       
      }
      return(dists)
    }
    
    if (length(mcs_with_primary) > 0) {
      c_dists = sapply(mcs_with_primary, function(prim_mc) order_mgraph_mcs_by_mc_dist(mgraph, prim_mc, rep(NA, ncol(mc@mc_fp))))
      dists = apply(c_dists, 1, min, na.rm=T)
      if (any(is.infinite(dists))) {
        message(glue("Found {sum(is.infinite(dists))} disconnected metacells, skiping..."))
        next
      }
      dists_cols = setNames(colorRampPalette(c('white', brewer.pal(n=9, 'YlOrRd')))(length(unique(dists))), as.character(sort(unique(dists))))
      
      # plot mc2d with dists as colour
      png(scfigs_fn(c_id, "mc2d_by_prim_dist", dir=scfigs_dir(c_id, "prim_mc_dist")), 600, 600)
      plot(mc2d@mc_x, mc2d@mc_y, pch=21, cex=3, main=glue("{foc_bc}: Distance from MC with primary"), xaxt='n', yaxt='n', xlab='', ylab='')
      fr = mgraph$mc1
      to = mgraph$mc2
      segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to])
      points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=3, bg=dists_cols[as.character(dists)])
      text(mc2d@mc_x, mc2d@mc_y, names(mc2d@mc_x), cex=0.7)
      legend("topleft", legend=names(dists_cols), pch=21, pt.bg=dists_cols, bty='n', cex=0.8)
      dev.off()

      png(scfigs_fn(c_id, "mc2d_by_prim_dist_noMCid", dir=scfigs_dir(c_id, "prim_mc_dist")), 600, 600)
      plot(mc2d@mc_x, mc2d@mc_y, pch=21, cex=3, main=glue("{foc_bc}: Distance from MC with primary"), xaxt='n', yaxt='n', xlab='', ylab='')
      fr = mgraph$mc1
      to = mgraph$mc2
      segments(mc2d@mc_x[fr], mc2d@mc_y[fr], mc2d@mc_x[to], mc2d@mc_y[to])
      points(mc2d@mc_x, mc2d@mc_y, pch=21, cex=3, bg=dists_cols[as.character(dists)])
      legend("topleft", legend=names(dists_cols), pch=21, pt.bg=dists_cols, bty='n', cex=0.8)
      dev.off()
      
      # Modules correlated with prim dists
      lfp = log2(mc@mc_fp)
      gmax = apply(lfp, 1, max)
      foc_gs = names(which(gmax >= min_lfp_max))
      
      gs_d_cor = cor(t(lfp[foc_gs, ]), dists)[,1]
      top_cor_gs = names(tail(sort(gs_d_cor), n=disp_n_d_cor_gs))
      top_acor_gs = names(head(sort(gs_d_cor), n=disp_n_d_cor_gs))
      
      gs_d_cor_top_df = reshape2::melt(lfp[c(top_cor_gs, top_acor_gs), ], varnames=c('gene', 'mc_id'), value.name = 'lfp') %>%
        mutate(prim_dist = dists[mc_id],
               d_cor = ifelse(gs_d_cor[as.character(gene)] > 0, 'Pos', 'Neg'),
               gene_ord = factor(gene, levels=c(rev(top_cor_gs), top_acor_gs)))
      
      p = ggplot(gs_d_cor_top_df, aes(x=prim_dist, y=lfp, fill=d_cor)) + geom_jitter(position=position_jitter(0.2), pch=21, size=0.8, show.legend=F) + labs(title=glue("{foc_bc}: top/bottom genes cor to prim dists"), x="", y="Enr (log2)") + stat_summary(fun = mean, fun.min = mean, fun.max = mean, color="black", geom = "crossbar", show.legend=F) + facet_wrap(~gene_ord, scales='free_y', nrow=2) + scale_fill_manual(values=c(Pos='yellow', Neg='blue'))
      ggsave(scfigs_fn(c_id, "top_genes_cor_with_dists", dir=scfigs_dir(c_id, "prim_mc_dist")), p, width=disp_n_d_cor_gs*2, height=3)
      
      # gene modules dynamics across dists 
      gg_cor = cor(t(lfp[foc_gs, ]))
      gg_hc = hclust(dist(gg_cor), method='ward.D2')
      gg_cls = cutree(gg_hc, h = mean(gg_hc$height) + n_hc_sd_for_cutree * sd(gg_hc$height))
      
      gg_cls_ann = data.frame(row.names = names(gg_cls), gm = paste0('GM', gg_cls))
      pal_nms = c('Set1', 'Set2', 'Set3', 'Paired', 'Dark2', 'Accent', 'Pastel1', 'Pastel2', 'Spectral', 'BrBG')
      gg_cls_cols = setNames(unique(unlist(sapply(pal_nms, function(nm) brewer.pal(n=brewer.pal.info[nm, 'maxcolors'], nm))))[1:length(unique(gg_cls))],
                             paste0('GM', 1:length(unique(gg_cls))))
      
      gg_cor_o = gg_cor[gg_hc$order, gg_hc$order]
      g_gaps = cumsum(rle(gg_cls[gg_hc$order])$lengths)
      
      png(scfigs_fn(c_id, "gene_modules_cor_hm", dir=scfigs_dir(c_id, "prim_mc_dist")), nrow(gg_cor) + 400, nrow(gg_cor) + 400)
      pheatmap(gg_cor_o, breaks=seq(-1, 1, len=101), cluster_rows=F, cluster_cols=F, cellwidth=1, cellheight=1, gaps_row=g_gaps, gaps_col=g_gaps, annotation_row=gg_cls_ann, annotation_col=gg_cls_ann, show_rownames=F, show_colnames=F, annotation_color=list(gm=gg_cls_cols))
      dev.off()
      
      tfs_cls = gg_cls[names(gg_cls) %in% tfs]
      gm_top_tfs = lapply(split(names(tfs_cls), tfs_cls), function(gs) paste0(names(head(sort(gmax[gs], decr=T), disp_n_d_cor_gs)), collapse="\n")) 
      
      non_tf_gg_cls = setdiff(names(gg_cls), names(tfs_cls))
      cls_strong_gs = tapply(gmax[non_tf_gg_cls], gg_cls[non_tf_gg_cls], function(v) paste0(names(tail(sort(v), disp_n_d_cor_gs)), collapse=', '))
      
      gm2gs = sapply(split(non_tf_gg_cls, gg_cls[non_tf_gg_cls]), function(gm_gs) paste0(names(sort(gmax[gm_gs], dec=T)), collapse=', '))
      gm2tfs = sapply(split(names(tfs_cls), tfs_cls), function(gm_gs) paste0(names(sort(gmax[gm_gs], dec=T)), collapse=', '))
      gm_df = data.frame(row.names = names(gm2gs), pass_filter=F, genes=gm2gs, TFs=NA, terms=NA)
      gm_df[names(gm2tfs), 'TFs'] = gm2tfs
      
      cls_mu_lfp_per_dist = tgs_matrix_tapply(tgs_matrix_tapply(lfp[names(gg_cls), ], dists, mean), gg_cls, mean)
      gm_mu_lfp = tgs_matrix_tapply(t(lfp[names(gg_cls), ]), gg_cls, mean)
      cls_gm_se_per_dist =  t(tgs_matrix_tapply(gm_mu_lfp, dists, sd, na.rm=T) / sqrt(tgs_matrix_tapply(gm_mu_lfp, dists, length)))
      cls_gm_se_per_dist[is.na(cls_gm_se_per_dist)] = 0
      
      if (is.null(glob_mu_gm_dist)) {
        glob_mu_gm_dist = cls_mu_lfp_per_dist
        rownames(glob_mu_gm_dist) = paste(foc_bc, rownames(glob_mu_gm_dist))
      } else {
        c_gm_mu_dists = cls_mu_lfp_per_dist
        rownames(c_gm_mu_dists) = paste(foc_bc, rownames(c_gm_mu_dists))
        if (ncol(glob_mu_gm_dist) <= ncol(cls_mu_lfp_per_dist)) {
          new_dists = setdiff(colnames(cls_mu_lfp_per_dist), colnames(glob_mu_gm_dist))
          glob_mu_gm_dist = cbind(glob_mu_gm_dist, matrix(NA, nrow=nrow(glob_mu_gm_dist), ncol=length(new_dists), dimnames=list(rownames(glob_mu_gm_dist), new_dists) ))
        } else {
          new_dists = setdiff(colnames(glob_mu_gm_dist), colnames(c_gm_mu_dists))  
          c_gm_mu_dists = cbind(c_gm_mu_dists, matrix(NA, nrow=nrow(c_gm_mu_dists), ncol=length(new_dists), dimnames=list(rownames(c_gm_mu_dists), new_dists) ))
        }
        
        glob_mu_gm_dist = rbind(glob_mu_gm_dist, c_gm_mu_dists)
      }
      
      gm_gsea_non_filt = do.call('rbind', lapply(rownames(cls_mu_lfp_per_dist), 
                                        function(gm) { res = gsea_wrapper(gene_list=sort(gmax[names(which(gg_cls == gm))], dec=T), scoreType='pos', pval=1, verbose=F); if(nrow(res) > 0) { res$gm = gm; res }}))
      gm_gsea_non_filt$leadingEdge = sapply(1:nrow(gm_gsea_non_filt), function(i) paste0(gm_gsea_non_filt[[i, 'leadingEdge']], collapse=", "))
      gm_gsea_top_term = group_by(gm_gsea_non_filt, gm) %>% slice_min(pval, n = 1)
      gm_gsea = filter(gm_gsea_non_filt, padj <= gsea_max_padj)
      if (nrow(gm_gsea) > 0) {
        write.csv(gm_gsea, scfigs_fn(c_id, glue("GM_gsea_enriched_terms_pvalLE{gsea_max_padj}"), dir=scfigs_dir(c_id, "prim_mc_dist"), ext='csv'), row.names = F)
        write.csv(gm_gsea_top_term, scfigs_fn(c_id, "GM_gsea_top_term_per_gm", dir=scfigs_dir(c_id, "prim_mc_dist"), ext='csv'), row.names = F)
        gm2terms = sapply(split(gm_gsea$pathway, gm_gsea$gm), paste0, collapse=', ')
        gm_df[names(gm2terms), 'terms'] = gm2terms
      }
      
      selected_gms = apply(apply(cls_mu_lfp_per_dist, 1, range), 2, diff) >= min_gm_mu_diff
      cls_mu_lfp_per_dist_f = cls_mu_lfp_per_dist[selected_gms, ]
      cls_gm_se_per_dist_f = cls_gm_se_per_dist[selected_gms, ]
      
      gm_df[rownames(cls_mu_lfp_per_dist_f), 'pass_filter'] = T
      write.csv(rownames_to_column(gm_df, var='gm'), scfigs_fn(c_id, "GM_genes_TFs_terms", dir=scfigs_dir(c_id, "prim_mc_dist"), ext='csv'), row.names = F)
      gm_df = rownames_to_column(gm_df, var='gm') %>%
        mutate(clone = foc_bc)
      glob_gm_df = rbind(glob_gm_df, gm_df)
      
      rownames(cls_mu_lfp_per_dist_f) = glue("{rownames(cls_mu_lfp_per_dist_f)}: {cls_strong_gs[rownames(cls_mu_lfp_per_dist_f)]}")
      cls_ord_by_mean_enr = order(as.vector(cls_mu_lfp_per_dist_f %*% 1:ncol(cls_mu_lfp_per_dist_f)))
      cmlpdf_o = cls_mu_lfp_per_dist_f[cls_ord_by_mean_enr, ]
      cmlpdf_se_o = cls_gm_se_per_dist_f[cls_ord_by_mean_enr, ]
      
      c_ord = order(apply(cmlpdf_o, 1, which.max) + 1e-3 * apply(cmlpdf_o, 1, max))
      cmlpdf_o = cmlpdf_o[c_ord, ]
      cmlpdf_se_o = cmlpdf_se_o[c_ord, ]
      
      tfs_cls = gg_cls[names(gg_cls) %in% tfs]
      gm_top_tfs = lapply(split(names(tfs_cls), tfs_cls), function(gs) paste0(names(head(sort(gmax[gs], decr=T), disp_n_d_cor_gs)), collapse="\n")) 
      
      zlim = 1; cs=20
      png(scfigs_fn(c_id, "gmod_dynamics_on_dists_hm", dir=scfigs_dir(c_id, "prim_mc_dist")), ncol(cls_mu_lfp_per_dist_f) * cs + 500, nrow(cls_mu_lfp_per_dist_f) * cs + 200)
      pheatmap(pmin(pmax(cmlpdf_o, -zlim), zlim), breaks=seq(-zlim, zlim, len=101), cluster_cols=F, cluster_rows=F, cellwidth=cs, cellheight=cs, main=glue("{foc_bc}: Gene modules dynamics across prim dists"))
      dev.off()    
      
      png(scfigs_fn(c_id, "gmod_dynamics_on_dists_b", dir=scfigs_dir(c_id, "prim_mc_dist")), ncol(cmlpdf_o) * 30 + 500, nrow(cmlpdf_o) * 75 + 450, res=150, pointsize=7)
      layout(matrix(c(seq(4, len=nrow(cmlpdf_o)), 1), ncol=1), heights = c(rep(1, nrow(cmlpdf_o)), 1.5))
      par(mar=c(1,8,2,8))
      
      epi = scdb_gset(epi_gset)@gene_set
      epi = epi[intersect(names(epi), mat@genes)]
      dist_tot_umis = as.vector(tapply(colSums(mat@mat[, names(mc@mc)]), dists[mc@mc], sum))
      dist_epi_tot_umis = tgs_matrix_tapply(tgs_matrix_tapply(mat@mat[names(epi), names(mc@mc)], dists[mc@mc], sum), epi, sum)
      dist_epi_f_umis = t(dist_epi_tot_umis) / dist_tot_umis
      
      y_expand = diff(range(dist_epi_f_umis)) / 10
      bp = barplot(t(dist_epi_f_umis[, 1]), ylim=c(min(dist_epi_f_umis) - y_expand, max(dist_epi_f_umis) + y_expand), col=NA, border=NA, ylab='%UMIs', xaxt='n', main='Normal epithelial signatures', xlab='')
      lapply(colnames(dist_epi_f_umis), function(nm) { lines(bp, dist_epi_f_umis[, nm], col=supmc_cols[nm]);  points(bp, dist_epi_f_umis[, nm], col=supmc_cols[nm], pch=19, cex=2) })
      par(xpd=T)
      lo = legend(max(bp) + 0.5, max(dist_epi_f_umis), legend=colnames(dist_epi_f_umis), fill=supmc_cols[colnames(dist_epi_f_umis)], bty='n')
      rect(par('usr')[1], min(dist_epi_f_umis) - y_expand, lo$rect$left + lo$rect$w, max(dist_epi_f_umis) + y_expand)
      par(xpd=F)

      par(mar=c(1,8,2,8))
      for (i in 1:nrow(cmlpdf_o)) { 
        c_gm = gsub(": .*", "", rownames(cmlpdf_o)[i])
        top_term = filter(gm_gsea, gm == c_gm) %>% arrange(padj) %>% head(n=1) %>% pull(pathway)
        gm_title = paste0('GM', c_gm, ifelse(length(top_term) > 0, paste0(": ", top_term), ""))
        
        c_ylim = range(rbind(cmlpdf_o[i, ] + 2 * cmlpdf_se_o[i, ], cmlpdf_o[i, ] - 2 * cmlpdf_se_o[i, ]), na.rm=T)
        c_ylim = range(c_ylim + c(-1, 1) * diff(c_ylim)/6)
        bp = barplot(cmlpdf_o[i,], col=NA, border=NA, xaxt='n', main=gm_title, xlab='', ylab='', ylim=c_ylim)
        lines(bp, cmlpdf_o[i, ])
        c_se = cmlpdf_se_o[i, ]
        f_se = !is.na(c_se)
        arrows(bp[f_se], cmlpdf_o[i, f_se] - 2 * c_se[f_se], bp[f_se], cmlpdf_o[i, f_se] + 2 * c_se[f_se], code=3, angle=90, length=0.075)
        points(bp, cmlpdf_o[i, ], pch=21, cex=2, bg=dists_cols)
        abline(h=0, lty=2)
        axis(4, at=mean(range(cmlpdf_o[i, ])), labels=gsub(".*: ", "", gsub(", ", "\n", rownames(cmlpdf_o)[i])), tick=F, las=2)
        if (c_gm %in% names(gm_top_tfs)) {
          axis(2, at=mean(range(cmlpdf_o[i, ])), labels=gm_top_tfs[[c_gm]], tick=F, las=2, font=2, line=2)
        }
      }
      
      dev.off()
      
      # All genes in a GM dynamics across distances
      gm_odir = paste0(scfigs_dir(c_id, "prim_mc_dist"), "/gm_genes_dynamics")
      dir.create(gm_odir, showWarnings = F)
      for (i in 1:nrow(gm_df)) {
        c_gm = gm_df$gm[i]
        c_gs = strsplit(gm_df$genes[i], split=", ")[[1]]
        c_tfs = strsplit(gm_df$TFs[i], split=", ")[[1]]
        
        c_all = c(c_gs, c_tfs)
        c_all = c_all[!is.na(c_all)]
        
        c_mu = tgs_matrix_tapply(lfp[c_all, ], dists, mean)
        
        n_gs = length(c_all)
        nr = min(n_gs, ceiling(n_gs / 4))
        nc = ceiling(n_gs / nr)
        
        
        png(glue("{gm_odir}/GM{c_gm}_dist_dynamics.png"), nc * 350, nr * 150, res=150, pointsize=7)
        layout(matrix(seq(1, nr * nc), ncol=nc, byrow=F))
        par(mar=c(1,4,0.5,1))
        for(c_g in c_all) { 
          plot(as.numeric(rownames(c_mu)), c_mu[, c_g], type='b', ylab='', xlab='', xaxt='n', col=ifelse(c_g %in% c_tfs, 'blue', 'black'))
          abline(h=0, lty=2) 
          mtext(c_g, side=2, line=2, cex=1.5)
        }
        dev.off()               
        
        dir.create(glue("{gm_odir}/GM{c_gm}"), showWarnings = F, recursive = T)
        for (c_g in c_all) {
          xx = setNames(dists_cols, paste0('d', names(dists_cols)))
          p = ggplot(data.frame(lfp=lfp[c_g, ], prim_dist=factor(paste0("d", dists), levels=paste0("d", seq(0, max(dists))))), aes(x=prim_dist, y=lfp, fill=prim_dist)) + geom_jitter(position=position_jitter(0.2), pch=21, size=0.8, show.legend=F) + labs(title=sprintf("%s: %s%s", foc_bc, c_g, ifelse(c_g %in% c_tfs, " (TF)", "")), x="", y="Enr (log2)") + stat_summary(fun = mean, fun.min = mean, fun.max = mean, color="black", geom = "crossbar", show.legend=F) + scale_fill_manual(values=xx)
          ggsave(glue("{gm_odir}/GM{c_gm}/{c_g}.png"), p, width=3.5, height=2)
        }
      }
      
      # TFs by dist
      tfs_cls = gg_cls[names(gg_cls) %in% tfs]
      
      tfs_dist_mu_lfp = tgs_matrix_tapply(lfp[names(tfs_cls), ], dists, mean)
      tfs_dist_mu_lfp_f = tfs_dist_mu_lfp[, apply(abs(tfs_dist_mu_lfp), 2, max) >= 0.5 & apply(apply(tfs_dist_mu_lfp, 2, range), 2, diff) >= 0.5]
      
      if (ncol(tfs_dist_mu_lfp_f) > 0) {
        nr = min(ncol(tfs_dist_mu_lfp_f), ceiling(ncol(tfs_dist_mu_lfp_f) / 4))
        nc = ceiling(ncol(tfs_dist_mu_lfp_f) / nr)
        
        tfs_hc = hclust(dist(t(tfs_dist_mu_lfp_f)), method='ward.D2')
        
        png(scfigs_fn(c_id, "tfs_dynamics", dir=scfigs_dir(c_id, "prim_mc_dist")), nc * 350, nr * 150, res=150, pointsize=7)
        layout(matrix(seq(1, nr * nc), ncol=nc, byrow=F))
        par(mar=c(1,4,0.5,1))
        for(c_tf in colnames(tfs_dist_mu_lfp_f)[tfs_hc$order]) { 
          plot(as.numeric(rownames(tfs_dist_mu_lfp_f)), tfs_dist_mu_lfp_f[, c_tf], type='b', ylab='', xlab='', xaxt='n', col='blue')
          abline(h=0, lty=2) 
          mtext(c_tf, side=2, line=2, cex=1.5)
        }
        
        dev.off()               
      }

    } else {
      message(glue("{foc_bc}: No metacells found with >= {min_cells_for_primary_mc[foc_bc]} cells (top mc has {max(n_primary_per_mc)} cells), skipping..."))
    }
  }
  
}

#
##
LN_clones_phenotypic_space = function(mat_id = "pdx_LN_v2_filt",
                                      mc_id = "pdx_LN_v2_filt",
                                      bc_info_csv = "LN_bc_info.csv",
                                      selected_models=NULL,
                                      foc_clones = c(Clone_1 = "STG139 P1 GFP_ACTGATATCGAGATCGAAACTGGTAAAACGT",
                                                     Clone_5 = "STG201 P1 GFP_ACTGAAATCGTGATGGAAAGCGGTAAAACAA",
                                                     Clone_274 = "NKI250 P1 GFP_ACTGTCATCGAGATCCAAAAGGGTAAAACTA",
                                                     Clone_271 = "IC07 P1 GFP_ACTGCCATCCTGATCCAAACCGGTTAAACAA",
                                                     Clone_101 = "AB040 P1 GFP_ACTGCCATCTGGATGCAAAAAGGTGAAACTG"),
                                      min_clone_size = 1, min_clone_size_for_sc_plots=10, scale_sig = T,
                                      pal_gset_id="pdx_LN_v2_filt_colByCluster_pdx_LN_v2",
                                      max_plots_per_row=6, plot_size=300, pt_cex=0.5)
{
  mat = scdb_mat(mat_id)
  bc_info = fread(bc_info_csv, sep=",")
  mc = scdb_mc(mc_id)

  odir = scfigs_dir(mc_id, "clones_ternary_plots")
  dir.create(odir, showWarnings = F)
  
  plot_legend("topleft", nms=names(md_cols$Subtype),   colors=md_cols$Subtype,   ofn=glue("{odir}/subtype_legend.png"), leg_title='Subtype')
  plot_legend("topleft", nms=names(md_cols$Replicate), colors=md_cols$Replicate, ofn=glue("{odir}/replicate_legend.png"), leg_title='Replicate')
  plot_legend("topleft", nms=names(md_cols$Model),     colors=md_cols$Model,     ofn=glue("{odir}/model_legend.png"), leg_title='Model')
  
  if (is.null(selected_models)) {
    selected_models = sort(unique(mat@cell_metadata$Model))
  }
  
  md = mat@cell_metadata[names(mc@mc), ] %>% 
    rownames_to_column(var='cellid') %>%
    left_join(bc_info) %>%
    mutate(subtype = model2subtype[Model]) %>%
    filter(!is.na(var_seq) & Model %in% selected_models)
  
  md_f = group_by(md, Model, var_seq) %>% 
    summarise(n=n()) %>%
    filter(n >= min_clone_size) %>%
    left_join(md)
  
  clone_types = group_by(md_f, var_seq) %>%
    summarise(reps = paste0(sort(unique(Replicate)), collapse=' ')) %>%
    mutate(in_P = grepl("P1", reps),
           in_S = grepl("S", reps),
           cln_type = ifelse(in_P, ifelse(in_S, "Propagating", "Transient"), ifelse(in_S, "Emerging", "P2")))
  
  md_f = left_join(md_f, dplyr::select(clone_types, var_seq, cln_type))
  
  pal_sig_gset = scdb_gset(pal_gset_id)@gene_set
  supmc2marks = lapply(split(names(pal_sig_gset), pal_sig_gset), intersect, y=mat@genes)
  
  foc_c = md$cellid
  model2cells = split(foc_c, mat@cell_metadata[foc_c, 'Model'])
  
  uc = colSums(mat@mat[, foc_c])
  y = sapply(supmc2marks, function(v)  { colSums(mat@mat[v, foc_c]) / uc })
  
  if (scale_sig) {
    y = scale(y, center=F)
  }
  
  y_n = y / rowSums(y)
  y_n = y_n[, c('ML', 'LP', 'Basal')]

  inner_plot_cells_per_clone_by_model = function(md_f, c_model, name, min_clone_size, group_cells_by, groups_colors) {
    c_bcs_tab = filter(md_f, Model == c_model) %>%
      group_by(var_seq) %>% 
      summarise(n_cells = n()) %>%
      filter(n_cells >= min_clone_size) %>%
      arrange(desc(n_cells))
    
    n_r = ceiling(nrow(c_bcs_tab) / max_plots_per_row)
    n_c = min(nrow(c_bcs_tab), max_plots_per_row)
    
    png(scfigs_fn(mc_id, glue("ternary_f_umis_per_epi_sig_{c_model}_{name}_min{min_clone_size}cells_{ifelse(scale_sig, 'scaledSig', 'rawSig')}"), dir=odir), n_c * plot_size, n_r * plot_size)
    par(mfrow=c(n_r, n_c), mar=c(2,2,2,2))
    
    for (c_bc in c_bcs_tab$var_seq) {
      c_bc_tab = filter(md_f, Model == c_model & var_seq == c_bc) %>%
        group_by(!!as.symbol(group_cells_by)) %>% 
        summarise(n_cells = n()) %>%
        arrange(desc(n_cells))
      
      plot.new()
      plot.window(xlim = c(-0.1, 1.1), ylim = c(0, sqrt(3)/2+0.1), asp=sqrt(3)/2)
      title(main=gsub("GFP_", "", c_bc), font=2)
      text(0, 0, colnames(y_n)[1], pos=2)
      text(1, 0, colnames(y_n)[2], pos=4)
      text(0.5, sqrt(3)/2, colnames(y_n)[3], pos=3)
      segments(c(0, 0, 0.5), c(0, 0, sqrt(3)/2), c(0.5, 1, 1), c(sqrt(3)/2, 0, 0))
      
      points(y_n[,2] + y_n[,3]/2, y_n[,3] * sqrt(3)/2, pch=19, cex=pt_cex, col='lightgray')
      f = model2cells[[c_model]]
      points(y_n[f,2] + y_n[f,3]/2, y_n[f,3] * sqrt(3)/2, pch=19, cex=pt_cex, col='darkgray')
      
      for (c_grp in c_bc_tab[, group_cells_by]) {
        c_cells = filter(md_f, Model == c_model & var_seq == c_bc & !!as.symbol(group_cells_by) == c_grp)
        c_y_n = y_n[c_cells$cellid, , drop=F]
        points(c_y_n[,2] + c_y_n[,3]/2, c_y_n[,3] * sqrt(3)/2, pch=19, col=groups_colors[c_grp], cex=pt_cex)
      } 
    }
    
    dev.off()
  }
    
  for (c_model in unique(md_f$Model)) {
    inner_plot_cells_per_clone_by_model(md_f, c_model, 'All', min_clone_size = min_clone_size_for_sc_plots, group_cells_by = 'Replicate', groups_colors = md_cols$Replicate) 
    inner_plot_cells_per_clone_by_model(filter(md_f, grepl("P", Replicate)), c_model, 'Primary', min_clone_size = min_clone_size_for_sc_plots, group_cells_by = 'Replicate', groups_colors = md_cols$Replicate[c('P1', 'P2')]) 
  }
  
  # Focus clones (Fig 3E)
  if (!is.null(foc_clones)) {
    png(scfigs_fn(mc_id, glue("ternary_f_umis_per_epi_sig_foc_clones_{ifelse(scale_sig, 'scaledSig', 'rawSig')}"), dir=odir), length(foc_clones) * plot_size, plot_size)
    par(mfrow=c(1, length(foc_clones)), mar=c(2,2,2,2))
    
    for (c_nm in names(foc_clones)) {
      c_model = strsplit(foc_clones[c_nm], split=' ')[[1]][1]
      c_reps = strsplit(foc_clones[c_nm], split=' ')[[1]][2]
      c_bc = strsplit(foc_clones[c_nm], split=' ')[[1]][3]
      
      c_bc_tab = filter(md_f, Model == c_model & grepl(c_reps, Replicate) & var_seq == c_bc) %>%
        group_by(Replicate) %>% 
        summarise(n_cells = n()) %>%
        arrange(desc(n_cells))
      
      plot.new()
      plot.window(xlim = c(-0.1, 1.1), ylim = c(0, sqrt(3)/2+0.1), asp=sqrt(3)/2)
      title(main=paste(c_model, c_reps, c_nm), font=2)
      text(0, 0, colnames(y_n)[1], pos=2)
      text(1, 0, colnames(y_n)[2], pos=4)
      text(0.5, sqrt(3)/2, colnames(y_n)[3], pos=3)
      segments(c(0, 0, 0.5), c(0, 0, sqrt(3)/2), c(0.5, 1, 1), c(sqrt(3)/2, 0, 0))
      
      points(y_n[,2] + y_n[,3]/2, y_n[,3] * sqrt(3)/2, pch=19, cex=pt_cex, col='lightgray')
      f = model2cells[[c_model]]
      points(y_n[f,2] + y_n[f,3]/2, y_n[f,3] * sqrt(3)/2, pch=19, cex=pt_cex, col='darkgray')
      
      c_cells = filter(md_f, Model == c_model & var_seq == c_bc)
      c_y_n = y_n[c_cells$cellid, , drop=F]
      points(c_y_n[,2] + c_y_n[,3]/2, c_y_n[,3] * sqrt(3)/2, pch=19, col=md_cols$Model[c_model], cex=pt_cex)
       
    }
    
    dev.off()
  }
  
  # mean position per clone
  inner_plot_mean_clone = function(md_f, name, min_clone_size, group_clones_by, groups_colors) {
    mu_cln = as.data.frame(y_n) %>% 
      rownames_to_column(var='cellid') %>% 
      inner_join(dplyr::select(md_f, var_seq, !!as.symbol(group_clones_by), cellid)) %>%
      group_by(var_seq, !!as.symbol(group_clones_by)) %>%
      summarise(ML=mean(ML), LP=mean(LP), Basal=mean(Basal), n=n()) %>%
      filter(n >= min_clone_size) %>%
      arrange(desc(n))
    
    png(scfigs_fn(mc_id, glue("ternary_f_umis_per_epi_sig_{name}_muCln_min{min_clone_size}cells_{ifelse(scale_sig, 'scaledSig', 'rawSig')}_byModel"), dir=odir), plot_size * 2, plot_size * 2)
    par(mar=c(2,2,2,2))
    
    plot.new()
    plot.window(xlim = c(-0.1, 1.1), ylim = c(0, sqrt(3)/2+0.1), asp=sqrt(3)/2)
    title(main=glue("Mean {name} clones"), font=2)
    text(0, 0, colnames(y_n)[1], pos=2)
    text(1, 0, colnames(y_n)[2], pos=4)
    text(0.5, sqrt(3)/2, colnames(y_n)[3], pos=3)
    segments(c(0, 0, 0.5), c(0, 0, sqrt(3)/2), c(0.5, 1, 1), c(sqrt(3)/2, 0, 0))
    
    points(mu_cln$LP + mu_cln$Basal/2, mu_cln$Basal * sqrt(3)/2, pch=21, cex=0.2 + log10(mu_cln$n), bg=groups_colors[pull(mu_cln, !!as.symbol(group_clones_by))])

    groups_colors = groups_colors[sort(unique(pull(mu_cln, !!as.symbol(group_clones_by))))]    
    legend("topleft", legend=names(groups_colors), pch=21, pt.bg=groups_colors, bty='n')
    cln_sz = seq(0, ceiling(log10(max(mu_cln$n))))
    legend("topright", legend=cln_sz, col='black', pch=19, pt.cex=0.2 + cln_sz, title = '#cells (log10)', bty='n', y.intersp=2, x.intersp=2)
    
    dev.off()
    
    invisible(mu_cln)
  }
  
  inner_plot_mean_clone(filter(md_f, grepl("P", Replicate)), "Primary_by_Model",   min_clone_size = min_clone_size, group_clones_by = "Model",   groups_colors = md_cols$Model) 
  inner_plot_mean_clone(filter(md_f, grepl("P", Replicate)), "Primary_by_Subtype", min_clone_size = min_clone_size, group_clones_by = "subtype", groups_colors = md_cols$Subtype) 
  
  mu_p = inner_plot_mean_clone(filter(md_f, grepl("P", Replicate) & cln_type %in% c('Propagating', 'Transient')), "Primary_by_CloneType", min_clone_size = min_clone_size, group_clones_by = "cln_type",   groups_colors = md_cols$CloneType) 
  mu_s = inner_plot_mean_clone(filter(md_f, grepl("S", Replicate) & cln_type %in% c('Propagating', 'Emerging')), "Secondary_by_CloneType", min_clone_size = min_clone_size, group_clones_by = "cln_type",   groups_colors = md_cols$CloneType) 
  
  # propagating clones: with arrows from primary to 2ary
  mu_sp = inner_join(mu_p, mu_s, by='var_seq')
  
  png(scfigs_fn(mc_id, glue("ternary_f_umis_per_epi_sig_Propagating_muCln_min{min_clone_size}cells_{ifelse(scale_sig, 'scaledSig', 'rawSig')}"), dir=odir), plot_size * 3, plot_size * 3, res=180, pointsize=6)
  par(mar=c(2,2,2,2))
  
  plot.new()
  plot.window(xlim = c(-0.1, 1.1), ylim = c(0, sqrt(3)/2+0.1), asp=sqrt(3)/2)
  title(main=glue("Mean propagating clones"), font=2)
  text(0, 0, colnames(y_n)[1], pos=2)
  text(1, 0, colnames(y_n)[2], pos=4)
  text(0.5, sqrt(3)/2, colnames(y_n)[3], pos=3)
  segments(c(0, 0, 0.5), c(0, 0, sqrt(3)/2), c(0.5, 1, 1), c(sqrt(3)/2, 0, 0))
  
  px = mu_sp$LP.x + mu_sp$Basal.x/2
  py = mu_sp$Basal.x * sqrt(3)/2
  sx = mu_sp$LP.y + mu_sp$Basal.y/2
  sy = mu_sp$Basal.y * sqrt(3)/2
  
  points(px, py, pch=21, cex=0.2 + log10(mu_sp$n.x)/2, bg=md_cols$Replicate['P1'])
  points(sx, sy, pch=21, cex=0.2 + log10(mu_sp$n.y)/2, bg=md_cols$Replicate['S1'])
  
  arrows(px, py, sx, sy, length=0.05)
  
  legend("topleft", legend=c('Primary', 'Secondary'), pch=21, pt.bg=c(md_cols$Replicate['P1'], md_cols$Replicate['S1']), bty='n')
  cln_sz = seq(0, ceiling(log10(max(c(mu_sp$n.y, mu_sp$n.x)))))
  legend("topright", legend=cln_sz, col='black', pch=19, pt.cex=0.2 + cln_sz/2, title = '#cells (log10)', bty='n', y.intersp=2, x.intersp=2)
  
  dev.off()
  
  png(scfigs_fn(mc_id, glue("ternary_f_umis_per_epi_sig_Propagating_muCln_min{min_clone_size}cells_{ifelse(scale_sig, 'scaledSig', 'rawSig')}_zoom"), dir=odir), plot_size * 3, plot_size * 3, res=180, pointsize=6)
  par(mar=c(2,2,2,2))
  
  plot.new()
  plot.window(xlim = c(-0.1, 0.1) + range(c(px, sx)), ylim = c(-0.1, 0.1) + range(c(py, sy)), asp=sqrt(3)/2)
  points(px, py, pch=21, cex=0.2 + log10(mu_sp$n.x)/2, bg=md_cols$Replicate['P1'])
  points(sx, sy, pch=21, cex=0.2 + log10(mu_sp$n.y)/2, bg=md_cols$Replicate['S1'])
  
  arrows(px, py, sx, sy, length=0.05)
  
  dev.off()
  
}
