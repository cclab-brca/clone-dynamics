require(R.utils)
require(data.table)
require(ggpubr)
require(rstatix)
require(tgutil)
require(glue)
library(dplyr)
library(fgsea)


####
# reload packages and re-init parameters and directories
rl = function(scdb_dir="scdb", scfigs_dir="figs", config_fn=NULL, force_init=T)
{
  if (!is.null(config_fn) && file.exists(config_fn)) {
    tgconfig::override_params(config_fn, package="metacell")
  }
  
  dir.create(scdb_dir, recursive = T, showWarnings = F)
  dir.create(scfigs_dir, recursive = T, showWarnings = F)
  if (!exists(".scdb_base") || force_init) {
    scdb_init(scdb_dir, force_reinit=force_init)
  }
  scfigs_init(scfigs_dir)
}

####
# Test if object exists
scdb_obj_exists = function(obj_type, obj_id) 
{
  if(!exists(".scdb")) {
    message("scdb not initialized")
    invisible(F)
  } else {
    invisible(file.exists(sprintf("%s/%s.%s.Rda", .scdb_base, obj_type, obj_id)))
  }
}


###
# return a named vector of colors to names (from the mc color_key table)
get_mc_col2group = function(mc, white_is_undet=T) {
  col2group = as.character(mc@color_key$group)
  names(col2group) = as.character(mc@color_key$color)
  col2group = col2group[ names(col2group) %in% unique(mc@colors)] 
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      col2group = c(col2group, c('white'='Undet'))
    }
  }
  col2group[unique(names(col2group))]
}

####
# return a named vector of colors to names (from the mc color_key table)
get_mc_group2col = function(mc, white_is_undet=T) {
  group2col = as.character(mc@color_key$color)
  names(group2col) = as.character(mc@color_key$group)
  group2col = group2col[ group2col %in% unique(mc@colors)] 
  if (white_is_undet) {
    if ('white' %in% mc@colors) {
      group2col = c(group2col, c('Undet'='white'))
    }
  }
  
  group2col[unique(names(group2col))]
}

###
# differential expression between 2 groups of metacell (if mcs1/2 supplied) or cells (if nms1/2 supplier)
# loo_grps - perform the comparison with leave-one-out and only keep genes that appear in all runs (and either enriched or depleted in all). loo_grps is a list of the loo groups mapping to cell names
diff_expr = function(mc, mat_ds, mcs1=NULL, mcs2=NULL, reg=5, min_max_umi=0, nms1=NULL, nms2=NULL, filter_outlier_genes=F, compare_top_mc_to_n_highest=3, max_top_to_n_highest_ratio=3, verbose=T, geo_mean=F, geo_mean_per_cell=F, enr_min_mean_umi=0.1, enr_min_f_pos=0.5, calculate_p_value=F, p_adjust_method='fdr', loo_grps=NULL)
{
  if (is.null(nms1)) {
    nms1 = names(mc@mc)[mc@mc %in% mcs1]
  }
  if (is.null(nms2)) {
    nms2 = names(mc@mc)[mc@mc %in% mcs2]
  }
  nms1 = intersect(colnames(mat_ds), nms1)
  nms2 = intersect(colnames(mat_ds), nms2)
  if (verbose) {
    message(sprintf("comparing %d vs %d cells", length(nms1), length(nms2)))
  }
  
  loo_de_run = function(c_loo_grp) {
    c_nms1 = nms1
    c_nms2 = nms2
    if (!is.null(c_loo_grp)) {
      c_nms1 = setdiff(nms1, loo_grps[[c_loo_grp]])
      c_nms2 = setdiff(nms2, loo_grps[[c_loo_grp]])
    }
    
    if (geo_mean) {
      df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), mu1=apply(mat_ds[, c_nms1], 1, function(y) {exp(mean(log(1+y)))-1}), mu2=apply(mat_ds[, c_nms2], 1, function(y) {exp(mean(log(1+y)))-1}), stringsAsFactors = F)
      df$tot1 = df$mu1 * length(c_nms1)
      df$tot2 = df$mu2 * length(c_nms2)
    } else {
      df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), tot1=Matrix::rowSums(mat_ds[, c_nms1]), tot2=Matrix::rowSums(mat_ds[, c_nms2]), stringsAsFactors = F)
    }
    df$mean1 = Matrix::rowMeans(mat_ds[, c_nms1])
    df$mean2 = Matrix::rowMeans(mat_ds[, c_nms2])
    df$f_pos1 = Matrix::rowMeans(mat_ds[, c_nms1] > 0)
    df$f_pos2 = Matrix::rowMeans(mat_ds[, c_nms2] > 0)
    
    norm_by = min(sum(df$tot1), sum(df$tot2))
    df$tot1 = df$tot1 / sum(df$tot1) * norm_by
    df$tot2 = df$tot2 / sum(df$tot2) * norm_by
    
    if (geo_mean && geo_mean_per_cell) {
      df$enr = log2( (df$mu1 + reg) / (df$mu2 + reg))
    } else {
      df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
    }
    
    
    # Filtering
    df = df[pmax(df$tot1, df$tot2) >= min_max_umi &
              ifelse(df$enr > 0, df$mean1, df$mean2) >= enr_min_mean_umi &
              ifelse(df$enr > 0, df$f_pos1, df$f_pos2) >= enr_min_f_pos, ]
    
    df = df[order(df$enr, decreasing=T), ]
    
    if (filter_outlier_genes) {
      fp = mc@mc_fp[intersect(rownames(mc@mc_fp), df$gene), ]
      if (!is.null(mcs1)) {
        fp = fp[, mcs1, drop=F]
      }
      gmax = apply(fp, 1, max)
      gnext = apply(fp, 1, function(v) head(tail(sort(v), n=compare_top_mc_to_n_highest), n=1) )
      df[rownames(fp), 'out_r'] =  gmax/gnext
      to_filt = !is.na(df$out_r) & df$out_r > max_top_to_n_highest_ratio & df$enr > 0
      if (sum(to_filt) > 0) {
        if (verbose) {
          message(sprintf("filtering %d outlier genes:", sum(to_filt)))
        }
        print(df[to_filt, ])
        df = df[!to_filt, ]
      }
    }
    return(df)
  }
  
  n_cores = tgconfig::get_param("mc_cores", "metacell")
  if (n_cores > 1) {
    doMC::registerDoMC(n_cores)
  }
  
  df = loo_de_run(NULL)
  if (!is.null(loo_grps)) {
    loo_retained_genes = do.call("rbind", alply(names(loo_grps), 1, loo_de_run, .parallel = n_cores > 1)) %>%
      group_by(gene) %>%
      summarise(n=n(), all_enr = all(enr > 0), all_dep = all(enr < 0)) %>%
      filter(n == length(loo_grps) & (all_enr | all_dep)) %>% 
      pull(gene)
    df = df[loo_retained_genes, ]
  }
  
  # Add p-value if requested (wilcox test per gene, correct p-value for multiple tests)
  if (calculate_p_value) {
    
    #all_pvals = setNames(unlist(plyr::alply(rownames(mat_ds), 1, function(g) { wilcox.test(mat_ds[g, c_nms1], mat_ds[g, c_nms2])$p.value }, .parallel = n_cores > 1)),	                       rownames(mat_ds))
    all_pvals = setNames(unlist(plyr::alply(df$gene, 1, function(g) { wilcox.test(mat_ds[g, nms1], mat_ds[g, nms2])$p.value }, .parallel = n_cores > 1)), df$gene)
    adj_pvals = p.adjust(all_pvals, method=p_adjust_method)
    
    df$pval = all_pvals[df$gene]
    df$pval_adjust = adj_pvals[df$gene]
  }
  
  invisible(df)	        
}

###
# General function to select feature genes, auto-compute parameres by default
select_feature_genes = function(gstat_id, gset_id=gstat_id, lateral_gset_id = NULL, T_vm = NULL, T_tot = NULL, T_top3 = NULL, verbose=T) 
{
  gs = scdb_gstat(gstat_id)
  if (is.null(T_vm)) {
    T_vm = quantile(gs$ds_vm_norm, 0.75)
  }
  if (is.null(T_tot)) {
    T_tot = round(median(gs$tot))
  }
  if (is.null(T_top3)) {
    T_top3 = round(median(gs$ds_top3))
  }
  if (verbose) {
    message(sprintf("marker selection filters for %s: \ntotal umis:\t%d (%.2f) > %.2f\nds top3:\t%d (%.2f) > %.2f\nds vm norm:\t%d (%.2f) > %.2f", gstat_id, sum(gs$tot > T_tot), mean(gs$tot > T_tot), T_tot, sum(gs$ds_top3 > T_top3), mean(gs$ds_top3 > T_top3), T_top3, sum(gs$ds_vm_norm > T_vm), mean(gs$ds_vm_norm > T_vm), T_vm))
  }	
  
  mcell_gset_filter_varmean(gstat_id, gset_id, T_vm=T_vm, force_new=T)
  mcell_gset_filter_cov(gstat_id, gset_id, T_tot=T_tot, T_top3=T_top3)
  mcell_plot_gstats(gstat_id, gset_id, max_vm=NULL)
  
  if (!is.null(lateral_gset_id)) {
    marker_gset = scdb_gset(gset_id)	
    lateral_gset = scdb_gset(lateral_gset_id)
    if (verbose){
      message(sprintf("removing %d lateral genes from markers, left with %d", length(intersect(names(marker_gset@gene_set), names(lateral_gset@gene_set))), length(setdiff(names(marker_gset@gene_set), names(lateral_gset@gene_set)))))
    }
    marker_gset = gset_new_restrict_gset(marker_gset, lateral_gset, inverse=T, "cgraph markers w/o lat genes")
    scdb_add_gset(gset_id, marker_gset)
  }
}

####
# Build blaclisted gene sets by gene names (mitochondrial, IG, ncRNA genes, RP-# genes and snoXXX genes)
meta_build_blacklist_gsets_by_gene_nms = function(all_id, ds_nm, gene_pref=NULL, mito_mode='minimal') 
{
  full_m = scdb_mat(all_id)
  
  # mitochondrial gene set
  mito_gset_id = sprintf("%s_mito", ds_nm)
  
  mt_cands = grep("^MT-", full_m@genes, v=T, perl=T, ignore.case = T)
  mt_genes = setNames(rep('neto MT', length(mt_cands)), mt_cands)
  if (mito_mode != 'minimal') {
    mt_cands = grep(paste0(paste0("^", gene_pref, c("MT-", "MTRN", "MTAT", "MTND", "MRP")), collapse="|"), full_m@genes, v=T, perl=T, ignore.case = T)
    mito = data.table::fread("/Users/eyal-l01/proj/mc_common/data/MitoCarta2_human.txt", header=T, sep="\t", stringsAsFactors=F)
    mito_genes = paste0(gene_pref, mito$Symbol)
    mt_both = intersect(mt_cands, mito_genes)
    mt_cands = setdiff(mt_cands, mt_both)
    mitocarta = setdiff(mito_genes, mt_both)
    mt_genes = c(rep('regexp MT', length(mt_cands)), rep('MitoCarta2', length(mitocarta)), rep('MitoCarta2 and regexp MT', length(mt_both)))
    names(mt_genes) = c(mt_cands, mitocarta, mt_both)
  }
  scdb_add_gset(mito_gset_id, gset_new_gset(mt_genes, 'mitochondrial genes'))
  
  # IG gene set
  ig_gset_id = sprintf("%s_ig", ds_nm)
  ig_nms = grep(paste0(paste0("^", gene_pref, c("IGK", "IGL", "IGJ", "IGH", "IGBP", "IGSF")), collapse="|"), full_m@genes, v=T, perl=T, ignore.case = T)
  ig_genes = rep("IG", length(ig_nms))
  names(ig_genes) = ig_nms
  scdb_add_gset(ig_gset_id, gset_new_gset(ig_genes, 'IG genes'))
  
  # ncRNA gene set
  ncrna_gset_id = sprintf("%s_ncrna", ds_nm)
  ncrna_nms = paste0(gene_pref, c('MALAT1', 'XIST', 'NEAT1', 'hsa-mir-6723'))
  ncrna_genes = rep("ncRNA", length(ncrna_nms))
  names(ncrna_genes) = ncrna_nms
  scdb_add_gset(ncrna_gset_id, gset_new_gset(ncrna_genes, "ncRNA genes"))
  
  # RP pseudo genes set
  ncrp_gset_id = sprintf("%s_ncrp", ds_nm)
  ncrp_nms = grep(paste0("^", gene_pref, "sRP[0-9]+-"), full_m@genes, v=T, perl=T, ignore.case = T)
  ncrp_genes = rep("ncRP", length(ncrp_nms))
  names(ncrp_genes) = ncrp_nms
  scdb_add_gset(ncrp_gset_id, gset_new_gset(ncrp_genes, "RP##- genes"))
  
  # sno Genes
  sno_gset_id = sprintf("%s_sno", ds_nm)
  sno_nms = grep(paste0("^", gene_pref, "SNOR[AD][0-9]+"), full_m@genes, v=T, perl=T, ignore.case = T)
  sno_genes = rep("sno", length(sno_nms))
  names(sno_genes) = sno_nms
  scdb_add_gset(sno_gset_id, gset_new_gset(sno_genes, "SNORA/D genes"))
  
}

#
# build clean mat: remove blaclisted genes (mitochondrial, ncRNA, RP[0-9], IG) and small cells
meta_build_blist_filtered_master_mat = function(all_id, filt_id, ds_nm, min_umis_post_gene_ignore = 500, min_umis_pre_gene_ignore=500, max_mito_f = 0.6, filt_mat_by_column=NULL, sample_field=NULL, gsets_to_filter=c("mito", "ig", "ncrna", "ncrp", "sno"))
{
  full_m = scdb_mat(all_id)
  stopifnot(is.null(filt_mat_by_column) || filt_mat_by_column %in% colnames(full_m@cell_metadata) || all(is.logical(full_m@cell_metadata[, filt_mat_by_column])))
  
  blist_gsets = sapply(gsets_to_filter, function(gnm) scdb_gset(sprintf("%s_%s", ds_nm, gnm)))
  blist_genes = unlist(lapply(blist_gsets, function(gs) names(gs@gene_set)))
  
  uc = Matrix::colSums(full_m@mat)
  
  mt_genes = intersect(names(blist_gsets[["mito"]]@gene_set), full_m@genes)
  mito_f = Matrix::colSums(full_m@mat[mt_genes, ]) / uc
  
  full_m@cell_metadata$f_mito = mito_f
  scdb_add_mat(all_id, full_m)
  
  # plot %mito vs log umis
  png(scfigs_fn(all_id, "fMito_vs_logUMIs"), 300, 300)
  valid_cols = ifelse(mito_f >= max_mito_f, 'red', 'black')
  if (!is.null(filt_mat_by_column)) {
    valid_cols = ifelse(full_m@cell_metadata[, filt_mat_by_column], 'black', 'red')
  }
  plot(log2(uc), mito_f, pch=19, cex=0.5, col=valid_cols , xlab="UMIs (log2)", ylab="%mito")
  dev.off()
  
  # plot %mito by sample (if field was supplied)
  if (!is.null(sample_field)) {
    png(scfigs_fn(all_id, "fMito_by_sample"), 500, 500)
    par(mar=c(16, 4, 1, 1))
    boxplot(mito_f ~ full_m@cell_metadata[, sample_field], las=2, col='navyblue', notch=T, pch=19, cex=0.5, xlab='', ylab="% mito")
    abline(h=max_mito_f, col='red', lty=2)
    dev.off()
    
  }
  # filter cells with low counts, large MT fraction, ignore MT genes, RP##- and mega-strong RNA genes
  mcell_mat_ignore_genes(filt_id, all_id, blist_genes)
  
  filt_mat = scdb_mat(filt_id)
  to_ignore = filt_mat@cells[mito_f >= max_mito_f | Matrix::colSums(filt_mat@mat) <= min_umis_post_gene_ignore | uc <= min_umis_pre_gene_ignore]
  if (!is.null(filt_mat_by_column)) {
    to_ignore = filt_mat@cells[!filt_mat@cell_metadata[, filt_mat_by_column]]
  }
  
  mcell_mat_ignore_cells(filt_id, filt_id, union(filt_mat@ignore_cells, to_ignore))
  
  # working on the filtered mat
  mcell_add_gene_stat(filt_id, filt_id)
}

####
# Wrapper to cluster genes based on the mat given the input anchor genes
select_gene_modules_by_anchor_genes = function(mat_id, gene_anchors, gset_nm, cor_thresh = 0.1, gene_anti=c(), sz_cor_thresh=0.7, nclusts=20, downsample_n=NA) 
{
  tab_fn = sprintf("%s/%s.txt", scfigs_dir(mat_id, "gmods_by_anchors"), gset_nm)
  message("mcell_mat_rpt_cor_anchors")
  
  mcell_mat_rpt_cor_anchors(mat_id=mat_id,
                            gene_anchors = gene_anchors,
                            cor_thresh = cor_thresh,
                            gene_anti = gene_anti,
                            tab_fn = tab_fn,
                            sz_cor_thresh=sz_cor_thresh)
  
  
  foc_gcor = read.table(tab_fn, sep="\t", h=T, stringsAsFactors=F, check.names=F)
  foc_neto = foc_gcor[,setdiff(colnames(foc_gcor),c("sz_cor","max","neg_max"))]
  foc_genes = if (is.null(ncol(foc_neto))) { setNames(rep(1, length(foc_neto)), rownames(foc_gcor)) } else { apply(foc_neto, 1, which.max) }
  gset = tgGeneSets(foc_genes, gset_nm)
  
  scdb_add_gset(gset_nm, gset)
  
  sub_mat_id = paste(mat_id, gset_nm, sep="_")
  
  mcell_mat_ignore_genes(sub_mat_id, mat_id, names(foc_genes), reverse=T)
  
  message("mcell_gset_split_by_dsmat")
  mcell_gset_split_by_dsmat(gset_nm, sub_mat_id, nclusts)
  mcell_plot_gset_cor_mats(gset_nm, sub_mat_id)
  
}

####
# Generate lateral gene sets (cluster by anchor genes and selects clusters containing the anchror genes) - cell cycle, IFN response and stress
meta_generate_lateral_gene_sets = function(filt_id, ds_nm, specie="human", anchors=NULL, nclusts=20, add_by_name=NULL, cor_thresh=0.1, gene_pref=NULL, selected_sets=NULL, rebuild=T, n_cores=1)
{
  if (is.null(anchors)) {
    anchors = switch(specie,
                     human = list(cell_cycle=c('MKI67', 'HIST1H1D', 'PCNA', 'SMC4', 'MCM3', 'TYMS', 'TOP2A', 'TUBB'), ifn=c('ISG15', 'OAS1', 'WARS', 'IFIT1'), stress=c("TXN", "HSP90AB1", "HSPA1A", "FOS", "JUN", "HIF1A"), hypoxia=c("VEGFA", "NDRG1", 'BNIP3')),
                     mouse = list(cell_cycle=c('Mki67', 'Hist1h1d', 'Pcna', 'Smc4', 'Mcm3', 'Tyms', 'Top2a'), ifn=c('Isg15', 'Wars', 'Ifit1'), stress=c("Hsp90ab1", "Hspa1a", "Fos", "Jun", "Hif1a")))
  }
  if (is.null(add_by_name)) {
    add_by_name = switch(specie, 
                         human = list(cell_cycle='^HIST|^CENP|^SMC[0-9]', ifn="^IFI"),
                         mouse = list(cell_cycle='^Hist|^Cenp|^Smc[0-9]', ifn="^Ifi"))
  }
  
  mat_f = scdb_mat(filt_id)
  
  if (is.null(selected_sets)) {
    selected_sets = names(anchors)
  }
  
  if (n_cores > 1 & length(selected_sets) > 1) {
    doMC::registerDoMC(cores=n_cores)
  }	
  
  inner_gset_by_anchor = function(gnm) {
    gs_id = paste(ds_nm, gnm, sep="_")
    gs_filt_id = paste(gs_id, "filt", sep="_")
    if (rebuild | !scdb_obj_exists("gset", gs_filt_id)) {
      message(sprintf("Processing %s...", gnm))
      
      curr_anchors = paste0(gene_pref, anchors[[gnm]])
      select_gene_modules_by_anchor_genes(filt_id, curr_anchors, gset_nm=gs_id, cor_thresh=cor_thresh, sz_cor_thresh=cor_thresh, nclusts=nclusts)
      
      gset = scdb_gset(gs_id)
      selected_clusts = unique(gset@gene_set[curr_anchors])
      for (cl in selected_clusts) {
        file.rename(sprintf("%s/%s_%s.gset_cors/%d.png", .scfigs_base, ds_nm, gnm, cl), sprintf("%s/%s_%s.gset_cors/filt_%d.png", .scfigs_base, ds_nm, gnm, cl))
      }
      mcell_gset_remove_clusts(gs_id, filt_clusts=selected_clusts, new_id = gs_filt_id, reverse=T)
      
      gset = scdb_gset(gs_filt_id)
      add_re = add_by_name[[gnm]]
      if (!is.null(add_re)) {
        gset = gset_add_genes(gset, setdiff(grep(add_re, mat_f@genes, v=T, perl=T), names(gset@gene_set)), max(gset@gene_set)+1)
      }
      scdb_add_gset(gs_filt_id, gset)		
    }
    gset = scdb_gset(gs_filt_id)
    names(gset@gene_set)
  }
  
  lat_genes = plyr::alply(selected_sets, 1, inner_gset_by_anchor, .parallel=n_cores > 1 & length(selected_sets) > 1)
  names(lat_genes) = selected_sets
  lat_genes_memb = unlist(lapply(seq_along(lat_genes), function(i) { v = lat_genes[[i]]; r = rep(i, length(v)); names(r) = v; r }))
  scdb_add_gset(paste0(ds_nm, "_lateral"), gset_new_gset(lat_genes_memb, sprintf('lateral: %s', paste0(names(anchors), collapse=", "))))
  
}

####
# Pipline generating mc from mat
meta_mat2mc = function(mat_id, cells=NULL, name="", lateral_gset_id=NULL, 
                       T_vm = NULL, T_tot = NULL, T_top3 = NULL, T_lfc = 3000, 
                       cgraph_knn = NULL, cgraph_downsamp=T, 
                       bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
                       mc_K=30, min_mc_size=30, mc_alpha=2, feat_gset_id=NULL,
                       split_mc_with_dbscan_and_filt_outliers=T, rebuild=T) 
{
  set.seed(42)
  
  # create submat if required
  if (!is.null(cells)) {
    stopifnot(nchar(name) > 0)
    new_id = paste(mat_id, name, sep="_")
    mcell_mat_ignore_cells(new_id, mat_id, cells, reverse=T)
    mcell_add_gene_stat(new_id, new_id)
  } else {
    new_id = mat_id
  }
  
  # Add gstats if it doesnt exist
  mcell_add_gene_stat(new_id, new_id, force=rebuild)
  
  # select genes to affect graph creation
  if (is.null(feat_gset_id)) {
    select_feature_genes(gstat_id = new_id, gset_id=new_id, lateral_gset_id = lateral_gset_id, T_vm = T_vm, T_tot = T_tot, T_top3 = T_top3) 
    feat_gset_id = new_id
  }
  
  # create cgraph
  if (is.null(cgraph_knn)) {
    mat = scdb_mat(new_id)	
    cgraph_knn = min(max(ceiling(mat@ncells / 200), 20), 250)
  }
  if (rebuild || !scdb_obj_exists('cgraph', new_id)) {
    message(sprintf("creating balanced Knn graph, K = %d", cgraph_knn))
    mcell_add_cgraph_from_mat_bknn(mat_id=new_id, gset_id=feat_gset_id, graph_id=new_id, K=cgraph_knn, dsamp=cgraph_downsamp)
  }
  message(new_id, " cgraph done")
  
  # bootstrap
  if (rebuild || !scdb_obj_exists('coclust', new_id)) {
    mcell_coclust_from_graph_resamp(new_id, new_id,
                                    min_mc_size=round(cgraph_knn/5),
                                    p_resamp=bootstrap_p_resamp,
                                    n_resamp=bootstrap_n_resamp)
  }
  message(new_id, " bootstrap done")
  
  # mc from coclust matrix
  if (rebuild || !scdb_obj_exists('mc', new_id)) {
    mcell_mc_from_coclust_balanced(new_id, new_id, new_id, K=mc_K, min_mc_size=min_mc_size, alpha=mc_alpha)
    
    # find and remove outliers
    # breaking down heterogenous metacells and removing outliers by extreme gene expression
    if (split_mc_with_dbscan_and_filt_outliers) {
      mcell_mc_split_filt(new_id, new_id, new_id, T_lfc=T_lfc, plot_mats=F)
      mcell_plot_outlier_heatmap(mc_id=new_id, mat_id = new_id, T_lfc  = T_lfc)
    }	
  }
}

###
# common mat to mc and plots pipline
common_pipe = function(ds_nm, build_mat_func, max_f_mit=0.6, min_umis_pre_gene_ignore=500, min_umis_post_gene_ignore=500, gsets_to_filter=c("mito", "ig", "ncrna", "ncrp", "sno"), filt_mat_by_column=NULL,
                       gset_cor_thresh=0.1, mc_cells=NULL, mc_name=NULL, 
                       T_vm = NULL, T_tot = NULL, T_top3 = NULL, T_lfc = 3000, 
                       cgraph_knn = NULL, cgraph_downsamp=T, 
                       bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
                       mc_K=30, min_mc_size=30, mc_alpha=2, 
                       metadata_fields = NULL, color_by_conf = F, col_by_cutree_k=0,
                       gene_pref=NULL, selected_sets=NULL,
                       sample_field=NULL, specie="human", annot_dfs=c(human="hg_Blueprint", mouse="mm_RNAseq"),
                       rebuild=T, out_base=ds_nm, filt_lateral_genes=T, fp_T_fold=2, min_confu_nmc=2, seed=42,
                       split_mc_with_dbscan_and_filt_outliers=T) 
{
  set.seed(seed)
  
  all_id <<- sprintf("%s_raw", ds_nm)
  filt_id <<- sprintf("%s_filt", ds_nm)
  lateral_gset_id = if (filt_lateral_genes) sprintf("%s_lateral", ds_nm) else NULL
  dir.create(ds_nm, showWarnings = F)
  
  if (!is.null(out_base)) {
    rl(scdb_dir=paste0(out_base, "/scrna_db"), scfigs_dir=paste0(out_base, "/figs"), config_fn=sprintf("config/%s.yaml", ds_nm))
  }
  
  # build raw mat
  if (!is.null(build_mat_func) & (rebuild | !(paste0("mat.", all_id) %in% scdb_ls("mat")))) {
    message(sprintf("%s: building raw mat...", ds_nm))
    do.call(build_mat_func, args=list(mat_id=all_id))
    scdb_init(.scdb_base, T)
  }
  
  if (!scdb_obj_exists("gstat", all_id)) {
    mcell_add_gene_stat(all_id, all_id)
  }
  min_umis_cutoff = mcell_plot_umis_per_cell(all_id)
  
  # filter cells and genes (similar to melanoma - mito, IG, ncRNA)
  if (rebuild | !(paste0("mat.", filt_id) %in% scdb_ls("mat"))) {
    message(sprintf("%s: filtering mat...", ds_nm))
    meta_build_blacklist_gsets_by_gene_nms(all_id, ds_nm, gene_pref=gene_pref)	
    meta_build_blist_filtered_master_mat(all_id, filt_id, ds_nm, max_mito_f=max_f_mit, min_umis_pre_gene_ignore = min_umis_pre_gene_ignore, min_umis_post_gene_ignore = min_umis_post_gene_ignore, filt_mat_by_column=filt_mat_by_column, sample_field=sample_field, gsets_to_filter = gsets_to_filter)
  }
  
  # build lateral gene sets
  if (filt_lateral_genes && (rebuild || !(paste0("gset.", lateral_gset_id) %in% scdb_ls("gset")))) {
    message(sprintf("%s: building lateral gene sets...", ds_nm))
    meta_generate_lateral_gene_sets(filt_id, ds_nm, specie=specie, cor_thresh = gset_cor_thresh, gene_pref=gene_pref, selected_sets=selected_sets)
  }
  
  # partition to metacells
  mc_id = sprintf("%s%s", filt_id, ifelse(is.null(mc_name), "", paste0("_", mc_name)))
  new_mat_id = mc_id
  if (rebuild | !(paste0("mc.", mc_id) %in% scdb_ls("mc"))) {
    message(sprintf("%s: partitioning to metacells...", ds_nm))
    meta_mat2mc(filt_id, lateral_gset_id = lateral_gset_id, cells=mc_cells, name=mc_name, 
                T_vm = T_vm, T_tot = T_tot, T_top3 = T_top3, T_lfc = T_lfc, 
                cgraph_knn = cgraph_knn, cgraph_downsamp = cgraph_downsamp, 
                bootstrap_n_resamp = bootstrap_n_resamp, bootstrap_p_resamp = bootstrap_p_resamp,
                mc_K = mc_K, min_mc_size = min_mc_size, mc_alpha = mc_alpha, split_mc_with_dbscan_and_filt_outliers= split_mc_with_dbscan_and_filt_outliers,
                rebuild=rebuild)
  }
  
}




#
#
# Marker genes per mc category
mc_get_marker_genes_per_mc_group = function(mc_id, mat_id, n_ds=NULL, 
                                            n_top_genes=Inf, 
                                            min_gene_max_lfp=0.5, min_avg_lfp_diff=0.25, 
                                            min_fold_enr=0.5, 
                                            enr_min_mean_umi=0.5, enr_min_f_pos = 0.8,  
                                            blacklist_gset=NULL, blacklist_re='^MT-|^RP[LS][0-9]', 
                                            restrict_to_groups=NULL, method="tot_group_umis",
                                            do_gsea=F, gsea_pval=0.05, msigdb_gsets=c('hallmark', 'GOBP'), n_top_terms_for_disp=3, 
                                            do_plot=F, plot_metric="mean_lfp", tfs_ifn = "/home/eyal-l01/proj/src/resources/Lambert2018_TF_names_v_1.01.txt", n_top_genes_for_disp=10, filter_multi_annot_markers=T, 
                                            run_name="", mat_ds=NULL)
{
  tfs = fread(tfs_ifn, header=F)$V1
  
  mc = scdb_mc(mc_id)
  
  col2grp = get_mc_col2group(mc)
  grp2col = get_mc_group2col(mc)
  
  if (is.null(restrict_to_groups)) {
    restrict_to_groups = names(grp2col)
  }  
  
  mat = scdb_mat(mat_id)
  
  blacklist_genes = NULL
  if (!is.null(blacklist_re)) {
    blacklist_genes = grep(blacklist_re, mat@genes, perl=T, v=T)  
  }
  
  if (!is.null(blacklist_gset)) {
    bl_gset_gs = names(scdb_gset(blacklist_gset)@gene_set)
    blacklist_genes = c(blacklist_genes, bl_gset_gs)
  }
  
  if (method == "tot_group_umis") {
    mat_f = scm_ignore_cells(mat, names(mc@mc), reverse=T)
    if (is.null(n_ds)) {
      n_ds = scm_which_downsamp_n(mat_f)
    }
    message(sprintf("Downsampling %d UMIs from mat...", n_ds))
    if (is.null(mat_ds)) {
      mat_ds = scm_downsamp(mat_f@mat, n_ds)
    }
    mat_ds = mat_ds[setdiff(rownames(mat_ds), blacklist_genes), ]
  }
  
  lfp = log2(mc@mc_fp)
  lfp = lfp[setdiff(rownames(lfp), blacklist_genes), ]
  lfp = lfp[apply(lfp, 1, max) >= min_gene_max_lfp, ]
  
  annot_to_mcs = split(1:ncol(mc@mc_fp), col2grp[mc@colors])
  stopifnot(all(restrict_to_groups %in% names(annot_to_mcs)))
  
  restrict_to_mcs = unlist(annot_to_mcs[restrict_to_groups])
  
  df = NULL
  gsea_df = NULL
  annot2top_terms = as.list(setNames(rep('None', length(restrict_to_groups)), restrict_to_groups))
  annot2top_genes = as.list(setNames(rep('None', length(restrict_to_groups)), restrict_to_groups))
  annot2top_tfs   = as.list(setNames(rep('None', length(restrict_to_groups)), restrict_to_groups))
  
  for (annot in restrict_to_groups) {
    message(sprintf("%s vs rest...", annot))
    foc_mcs = annot_to_mcs[[annot]]
    rest_mcs = setdiff(restrict_to_mcs, foc_mcs)
    
    if (method == "tot_group_umis") {
      x = diff_expr(mc, mat_ds, mcs1=foc_mcs, mcs2=rest_mcs, enr_min_mean_umi = enr_min_mean_umi, enr_min_f_pos = enr_min_f_pos) %>% filter(enr >= min_fold_enr)
      x$annot = annot
      x = head(x, n=min(nrow(x), n_top_genes))
    } else if (method == "mean_mcs_lfp") {
      x_v = sort(rowMeans(lfp[, foc_mcs, drop=F]) - rowMeans(lfp[, rest_mcs, drop=F]), decreasing = T)
      x_v = x_v[x_v >= min_avg_lfp_diff]
      if (length(x_v) > 0) {
        x = head(data.frame(annot=annot, gene=names(x_v), mean_mcs_diff=x_v), min(nrow(x), n_top_genes))
      }
    } else {
      stop(sprintf("Unknown method to DE (%s), expecting tot_group_umis or mean_mcs_lfp", method))
    }
    
    if (nrow(x) > 0) {
      df = rbind(df, x)
    } 
    else {
      warning(sprintf("No enriched genes found for %s...", annot))
    }
  }
  
  full_df = df
  
  if (filter_multi_annot_markers) {
    unique_markers = names(which(table(df$gene) == 1))
    df = filter(df, gene %in% unique_markers)
  }
  
  for (c_annot in restrict_to_groups) {
    x = filter(df, annot == c_annot)
    if (do_gsea && nrow(x) > 0) {
      gsea_res = gsea_wrapper(gene_list=sort(setNames(x$enr, x$gene), dec=T), msigdb_gsets=msigdb_gsets, scoreType='pos', pval=gsea_pval, verbose=F)
      
      if(nrow(gsea_res) > 0) {
        gsea_df = rbind(gsea_df, 
                        mutate(gsea_res,
                               leadingEdge = paste0(leadingEdge[[1]], collapse=", "),
                               annot = annot))
        annot2top_terms[[c_annot]] = gsea_res$pathway[seq(1, min(nrow(gsea_res), n_top_terms_for_disp))]
        
      }        
    }
    x_nonTF = filter(x, ! gene %in% tfs)
    x_TF = filter(x, gene %in% tfs)
    
    if (nrow(x_nonTF) > 0) {     
      annot2top_genes[[c_annot]] = head(x_nonTF$gene, min(nrow(x_nonTF), n_top_genes_for_disp))
    }
    if (nrow(x_TF) > 0) {     
      annot2top_tfs[[c_annot]] = head(x_TF$gene, min(nrow(x_TF), n_top_genes_for_disp))
    }
  }
  
  if (!is.null(gsea_df) && nrow(gsea_df) > 0) {
    write.csv(gsea_df, scfigs_fn(mc_id, glue("{run_name}_fgsea_pvalLE{gsea_pval}"), ext='csv'), row.names = F)
  }
  
  if (do_plot) {
    
    if (plot_metric == "mean_lfp") {
      annot_mu = tgs_matrix_tapply(t(lfp[df$gene, restrict_to_mcs]), df$annot, mean)
    } else if (plot_metric == 'f_umis') {
      annot_mu_c = t(tgs_matrix_tapply(t(mat@mat[df$gene, names(mc@mc)]), df$annot, sum)) / colSums(mat@mat[, names(mc@mc)])
      annot_mu = t(tgs_matrix_tapply(t(annot_mu_c), mc@mc, mean))
      annot_mu = annot_mu[, restrict_to_mcs]
    }
    write.csv(df, scfigs_fn(mc_id, glue("{run_name}_unique_sig_genes"), ext='csv'), row.names = F)
    write.csv(annot_mu, scfigs_fn(mc_id, glue("{run_name}_{plot_metric}_markers_sig"), ext='csv'), row.names = T, col.names=T)
    
    nms = sort(rownames(annot_mu))
    png(scfigs_fn(mc_id, glue("{run_name}_{plot_metric}_markers_sig_plts")), 250 * length(nms), 250 * (length(nms) + 1), res=180, pointsize=7)
    par(mfrow=c(length(nms) + 1, length(nms)), mar=c(4,4,2,2))
    for (i in 1:length(nms)) {
      for (j in 1:length(nms)) {
        if (i == j) {
          plot.new()
          plot.window(0:1, 0:1)
          box()
          legend("topleft", nms[i], pch=21, pt.bg = grp2col[nms[i]], bty='n', cex=0.8, pt.cex = 2)
          if (annot2top_terms[nms[i]] != 'None') {
            text(0, 0.5, paste0(annot2top_terms[nms[i]], collapse="\n"), pos=4, cex=0.5)
          }
        } else {
          plot(annot_mu[nms[i], ],  annot_mu[nms[j], ], pch=19, col=mc@colors[restrict_to_mcs], cex=1, xlab=nms[i], ylab=nms[j])
          #text(annot_lfp[nms[i], ], annot_lfp[nms[j], ], restrict_to_mcs, cex=0.7)
        }
      }
    }
    
    for (nm in nms) {
      plot.new()
      plot.window(0:1, 0:1)
      box()
      title(glue("{nm} ({length(annot2top_genes[[nm]])})"))
      legend("topleft", legend=annot2top_tfs[[nm]], bty='n', title='TFs', title.font=2)
      legend("topright", legend=annot2top_genes[[nm]], bty='n', title='non-TFs', title.font=2)
    }
    #mtext(run_name, 2, outer = T, cex = 2, font = 2)
    dev.off()
    
  }
  
  invisible(full_df)
  
}

# Adapted from: https://bioinformaticsbreakdown.com/how-to-gsea/
# Adjustments: removing nperm to use fgseaMultilevel
# If gene_list is not fold change but positive enrichment, set scoreType='pos'
####
gsea_wrapper = function(gene_list, msigdb_gsets=c('hallmark', 'GOBP'), msigdb_idir="/home/eyal-l01/proj/src/resources/MSigDB/msigdb_v2023.1.Hs_files_to_download_locally/msigdb_v2023.1.Hs_GMTs/", 
                        pval=0.05, minSize=15, maxSize=400, scoreType='std', plot_ofn=NULL, collapse_pathways=T, verbose=T) {
  set.seed(42)
 
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  nm2gmt = c(hallmark = "h.all.v2023.1.Hs.symbols.gmt",
             GO = "c5.go.v2023.1.Hs.symbols.gmt",
             GOBP = "c5.go.bp.v2023.1.Hs.symbols.gmt",
             reactome = "c2.cp.reactome.v2023.1.Hs.symbols.gmt",
             cell_types = "c8.all.v2023.1.Hs.symbols.gmt",
             all = "msigdb.v2023.1.Hs.symbols.gmt"
  )
  stopifnot(all(msigdb_gsets %in% names(nm2gmt)))
  
  myGO = do.call('c', lapply(msigdb_gsets, function(nm) fgsea::gmtPathways(paste0(msigdb_idir, nm2gmt[nm]))))
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize = minSize, ## minimum gene set size
                        maxSize = maxSize, ## maximum gene set size
                        #nperm = nperm,
                        scoreType = scoreType) %>% 
    as.data.frame()
  
  best_padj = min(fgRes$padj, na.rm=T)
  
  fgRes <- dplyr::filter(fgRes, padj < !!pval) %>% 
    arrange(desc(NES))
  
  if (verbose) { 
    message(glue("Number of signficant gene sets = {nrow(fgRes)} (best padj pre filter: {best_padj})"))
  }
  
  if (nrow(fgRes) > 0) {
    
    if (collapse_pathways) {
      if (verbose) { 
        message("Collapsing Pathways -----")
      }
      concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                          pathways = myGO,
                                          stats = gene_list)
      fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
      if (verbose) { 
        message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
      }
    }
    fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
    
    if (!is.null(plot_ofn)) {
      filtRes = rbind(head(fgRes, n = 10),
                      tail(fgRes, n = 10 ))
      
      total_up = sum(fgRes$Enrichment == "Up-regulated")
      total_down = sum(fgRes$Enrichment == "Down-regulated")
      header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
      
      colos = setNames(c("firebrick2", "dodgerblue2"),
                       c("Up-regulated", "Down-regulated"))
      
      p = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
        geom_point( aes(fill = Enrichment, size = size), shape=21) +
        scale_fill_manual(values = colos ) +
        scale_size_continuous(range = c(2,10)) +
        geom_hline(yintercept = 0) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title=header)
      ggsave(plot_ofn, p, width=10, height=nrow(filtRes) + 2)
    }
    
  } 
  
  return(fgRes)
}
