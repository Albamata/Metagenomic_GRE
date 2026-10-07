# ==============================================================================
#  Exploring the Bacteria Extracellular Vesicles fraction of 
#  the Guadalquivir River estuary (GRE)
#  Master's thesis (TFM)
#  Author: Alba Mata Gonzalez
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, readxl, phyloseq, vegan, patchwork, iNEXT)

options(dplyr.summarise.inform = FALSE, width = 130)

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
OPTIONAL_PKGS <- c("microbiome", "eulerr", "clusterProfiler", "KEGGREST")
missing_opt   <- OPTIONAL_PKGS[!vapply(OPTIONAL_PKGS, has_pkg, logical(1))]
if (length(missing_opt))
  message("Optional packages not installed (their blocks are skipped): ",
          paste(missing_opt, collapse = ", "))

# ---- 0.1 Configuration -------------------------------------------------------
CFG <- list(dir_out = "outputs",
            save_figs = TRUE,
            run_new = TRUE,      # run the additional checks tagged [NEW]
            seed = 123,
            min_pairs = 1,       # Fig. 4D: pairs in which a pathway must occur in both
  # fractions. The source script used 1; 3 is safer.
  font  = "sans",    # [FIX] "Helvetica" broke family names on export
  # ("Bur kholder iaceae" in Fig. 2)
  harmonise_colours = FALSE,   # TRUE = one zone and one fraction palette everywhere
  # --- inputs: mOTUs (PART A)
  in_motus_counts = "MetaQVIR_motus_filtered.txt",
  in_motus_tax = "MetaQVIR_motus_GTDB_taxonomy.tsv",
  in_motus_meta = "analisis_metada.txt",
  # --- inputs: SqueezeMeta tables (PARTS B and C)
  in_kegg_tpm = "../GUADALQUIVIR/resultado/kegg_tpm.tsv",
  in_kegg_tax = "../GUADALQUIVIR/resultado/kegg_tax.tsv",
  in_env = "../metadata_2024_11_11.xlsx",
  # --- inputs: KEGG pathway tables saved from KEGGREST (PART C)
  in_ko_path = "ko_path_tbl.rds",
  in_path_names = "pathway_names_tbl.rds",
  in_ko_info = "ko_info_tbl.rds",
  # --- inputs: NCycDB (PART D)
  in_ncyc = "Resultado_FC_only_PN_plus_FV.xlsx",
  in_ncyc_sheet = "Fusion_FC_PN_FC_FV_PN",
  in_ncyc_sets = "Ncyccompleteness.xlsx")

if (exists("CFG_OVERRIDE")) CFG <- modifyList(CFG, CFG_OVERRIDE)

dir.create(CFG$dir_out, showWarnings = FALSE, recursive = TRUE)
required_inputs <- unlist(CFG[c("in_motus_counts", "in_motus_tax", "in_motus_meta",
                                "in_kegg_tpm", "in_kegg_tax", "in_env",
                                "in_ncyc", "in_ncyc_sets")])
if (any(!file.exists(required_inputs)))
  stop("Input file(s) not found:\n  ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = "\n  "))

# ---- 0.2 Sampling design: single source of truth -----------------------------
# Two naming systems coexist in the raw files:
#   motus_id  A1 ... D3 / V_A1 ... V_D2     (mOTUs tables)
#   sqm_id    Z4_2 ... Z1_15 / Z4_V2 ...    (SqueezeMeta tables; Z4 = marine end)

ZONE_LEVELS <- c("Coastal Waters", "Inlet", "Transition", "Inner")
GROUP_LEVELS <- c("Marine-influenced", "Low-salinity")   # salinity > 10 vs < 1
FRAC_LEVELS <- c("Cellular", "Vesicular")

DESIGN <- tibble(
  motus_id = c("A1", "B1", "B2", "B3", "C1", "C2", "C3", "D1", "D2", "D3",
               "V_A1", "V_B2", "V_C2", "V_D2"),
  sqm_id = c("Z4_2", "Z3_5", "Z3_6", "Z3_7", "Z2_9", "Z2_10", "Z2_11",
               "Z1_13", "Z1_14", "Z1_15", "Z4_V2", "Z3_V6", "Z2_V10", "Z1_V14"),
  station = c(2, 5, 6, 7, 9, 10, 11, 13, 14, 15, 2, 6, 10, 14),
  fraction = factor(rep(FRAC_LEVELS, c(10, 4)), levels = FRAC_LEVELS),
  zone = factor(c("Coastal Waters", rep(c("Inlet", "Transition", "Inner"), each = 3),
                  ZONE_LEVELS), levels = ZONE_LEVELS)) %>%
  mutate(label = if_else(fraction == "Vesicular", paste0("V", station), as.character(station)),
         group = factor(if_else(zone %in% c("Coastal Waters", "Inlet"),
                                GROUP_LEVELS[1], GROUP_LEVELS[2]), levels = GROUP_LEVELS),
         x_idx = match(station, sort(unique(station))))

SAMPLE_ORDER <- DESIGN$label
CF_LAB <- DESIGN$label[DESIGN$fraction == "Cellular"]
VF_LAB <- DESIGN$label[DESIGN$fraction == "Vesicular"]

# Design keyed by station label (mOTUs part) and by SqueezeMeta id (KEGG part)
DESIGN_LAB <- DESIGN %>%
  dplyr::select(Label = label, Fraction = fraction, Zone = zone, Group = group,
                StationNum = station, x_idx)
DESIGN_SQM <- DESIGN %>%
  dplyr::select(Station = sqm_id, Label = label, Fraction = fraction, Zone = zone,
                Group = group, StationNum = station, x_idx)

# The four CF-VF pairs, marine -> inner
PAIRS <- DESIGN %>%
  dplyr::select(station, zone, fraction, sqm_id, label) %>%
  pivot_wider(names_from = fraction, values_from = c(sqm_id, label)) %>%
  filter(!is.na(sqm_id_Vesicular)) %>%
  transmute(station, zone, cf = sqm_id_Cellular, vf = sqm_id_Vesicular,
            cf_lab = label_Cellular, vf_lab = label_Vesicular) %>%
  arrange(zone)

# ---- 0.3 Palettes, labels and themes -----------------------------------------
PAL <- list(   # colours exactly as in the current figures
  zone_col = c("Coastal Waters" = "#BB5FD0", "Inlet" = "#6D98BE",     # Fig. 2D
               "Transition" = "#7AB59C", "Inner" = "#E08866"),        # Fig. 3  
  frac_tax = c(Cellular = "#C0392B", Vesicular = "#27AE60"),          # Fig. 2C
  frac_col  = c(Cellular = "#2C7FB8", Vesicular = "#DE7065"), # Fig. 4 and Fig. 5
)

FRAC_SHAPES <- c(Cellular = 16, Vesicular = 17)
BASE_COLS   <- c("#4472C4", "#ED7D31", "#E74C3C", "#70AD47", "#A9D18E", "#8064A2",
                 "#FFC000", "#26A69A", "#C0785A", "#F4A8A0",
                 "#1F77B4", "#B5651D", "#9E480E", "#43682B", "#636363", "#997300",
                 "#255E91", "#D4A6C8", "#86BCB6", "#B07AA1")
COL_OTHERS <- "#D9D9D9"; COL_UNCLASS <- "#7F7F7F"

# Axis and legend titles kept in one place. The thesis rubric bans abbreviations
# in figure legends, so the labels below avoid them where possible.
LAB <- list(
  stations  = "Stations",
  abs_tax = "Absolute abundance (marker-gene counts)",
  rel = "Relative abundance (%)",
  hill_tax = "Alpha diversity (Hill number, order 1)",
  abs_ko = "Reads assigned to KEGG orthologues (millions)",
  kegg_cat = "KEGG category",
  ko_rich = "KEGG orthologue richness",
  ko_hill = "Hill number, order 1",
  ko_n = "Number of KEGG orthologues",
  orf_n = "Open reading frames with KEGG orthologue",
  ratio = expression(log[2] * "(vesicular / cellular), paired mean"),
  n_rich = "Nitrogen-cycle gene families",
  n_conc = "Concentration relative to maximum (%)",
  n_abund = "Abundance relative to maximum (%)",
  n_compl = "Gene families detected (%)")

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text = element_text(family = CFG$font, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size + 1, colour = "black"),
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(2, "pt"),
      legend.title  = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = unit(8, "pt"),
      legend.key.width = unit(10, "pt"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(4, 4, 4, 4, "pt"))}

theme_tfm <- function(base_size = 9) {   
  theme_classic(base_size = base_size) +
    theme(text = element_text(family = CFG$font),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.key.size = unit(0.4, "cm"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.ticks = element_line(colour = "black"))}

TAG_THEME <- theme(plot.tag = element_text(face = "bold", size = 12))

# ---- 0.4 Helpers --------------------------------------------------------------
MS <- list()
ms <- function(key, value, digits = 4) {
  if (is.numeric(value)) value <- format(signif(value, digits), big.mark = ",", trim = TRUE)
  value <- paste(value, collapse = " | ")
  MS[[key]] <<- value
  cat(sprintf("[MS] %-58s %s\n", key, value))
  invisible(value)}

banner <- function(txt) cat("\n", strrep("=", 78), "\n ", txt, "\n", strrep("=", 78), "\n", sep = "")

save_fig <- function(plot, name, width, height) {
  if (!isTRUE(CFG$save_figs)) return(invisible(NULL))
  # The rubric asks for high-quality JPG files; the PDF is kept for the journal.
  for (ext in c("jpg", "pdf")) {
    f <- file.path(CFG$dir_out, paste0(name, ".", ext))
    tryCatch(ggsave(f, plot, width = width, height = height, units = "mm", dpi = 600),
             error = function(e) message("Could not save ", f, ": ", conditionMessage(e)))}}

save_tab <- function(x, name) {
  readr::write_tsv(as_tibble(x), file.path(CFG$dir_out, paste0(name, ".tsv")))}

# One-way PERMANOVA pseudo-F from a distance matrix
pseudo_F <- function(d, grp) {
  d <- as.matrix(d); grp <- as.factor(grp); n <- nrow(d); a <- nlevels(grp)
  ss_tot <- sum(d[lower.tri(d)]^2) / n
  ss_w <- sum(vapply(levels(grp), function(g) {
    i <- which(grp == g)
    if (length(i) < 2) return(0)
    dg <- d[i, i]
    sum(dg[lower.tri(dg)]^2) / length(i)
  }, numeric(1)))
  ((ss_tot - ss_w) / (a - 1)) / (ss_w / (n - a))}

# Exact permutation P for a two-group PERMANOVA by complete enumeration.
exact_p_two_groups <- function(d, grp) {
  grp <- as.factor(grp); stopifnot(nlevels(grp) == 2)
  n <- length(grp); k <- sum(grp == levels(grp)[1])
  f_obs <- pseudo_F(d, grp)
  combs <- combn(n, k)
  f_all <- apply(combs, 2, function(ix) { g <- rep("b", n); g[ix] <- "a"; pseudo_F(d, g) })
  c(F = f_obs, allocations = ncol(combs), P_exact = mean(f_all >= f_obs - 1e-12))}

# Good-Turing sample coverage 
good_turing <- function(x) {
  x <- x[x > 0]; n <- sum(x); f1 <- sum(x == 1); f2 <- sum(x == 2)
  if (f1 == 0) return(1)
  1 - (f1 / n) * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))}


# ==============================================================================
# PART A. TAXONOMY (mOTUs) -- Figure 2, Table 1, Figure S1, Section 3.2
# ==============================================================================

# ---- A1. Import and phyloseq object ------------------------------------------
otu <- read.table(CFG$in_motus_counts, sep = "\t", header = TRUE)
rownames(otu) <- otu[, 1]; otu[, 1] <- NULL
otu[] <- lapply(otu, function(x) suppressWarnings(as.numeric(as.character(x))))
colnames(otu) <- gsub("^X", "", colnames(otu))

motus_tax <- read.table(CFG$in_motus_tax, sep = "\t", header = TRUE)
rownames(motus_tax) <- motus_tax[, 1]; motus_tax[, 1] <- NULL
motus_tax[motus_tax == ""] <- NA

motus_meta <- read.table(CFG$in_motus_meta, sep = "\t", header = TRUE)
rownames(motus_meta) <- motus_meta[, 1]; motus_meta[, 1] <- NULL
motus_meta[motus_meta == ""] <- NA

ps_all <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE),
                   tax_table(as.matrix(motus_tax)),
                   sample_data(motus_meta))
ps_all <- prune_taxa(taxa_sums(ps_all) > 0, ps_all)

ps <- subset_taxa(ps_all, Kingdom %in% c("d__Bacteria", "d__Archaea"))
stopifnot(setequal(sample_names(ps), DESIGN$motus_id))

sd0 <- data.frame(sample_data(ps), check.names = FALSE)
if ("Fraction" %in% names(sd0)) {
  d_frac <- DESIGN$fraction[match(rownames(sd0), DESIGN$motus_id)]
  bad <- rownames(sd0)[(sd0$Fraction == "vesicles") != (d_frac == "Vesicular")]
  if (length(bad)) warning("DESIGN and metadata disagree on fraction for: ",
                           paste(bad, collapse = ", "))}

# Station numbers as sample names + design columns
sample_names(ps) <- DESIGN$label[match(sample_names(ps), DESIGN$motus_id)]
sd1 <- data.frame(sample_data(ps), check.names = FALSE)
d1  <- DESIGN_LAB[match(rownames(sd1), DESIGN_LAB$Label), ]
sd1$Fraction <- d1$Fraction; sd1$Zone <- d1$Zone; sd1$Group <- d1$Group
sd1$Label <- d1$Label
sample_data(ps) <- sample_data(sd1)

# ---- A2. Sequencing-depth diagnostics ------------------------------------------
depth <- tibble(Label = sample_names(ps), Depth = as.numeric(sample_sums(ps))) %>%
  left_join(DESIGN_LAB, by = "Label")
print(depth %>% dplyr::select(Label, Fraction, Zone, Depth))

# ---- A3. Composition: numbers quoted in the text -------------------------------
dom <- tax_glom(ps, taxrank = "Kingdom")
dom_counts <- as(otu_table(dom), "matrix")
rownames(dom_counts) <- as.vector(tax_table(dom)[, "Kingdom"])
ms("3.2 domain %, pooled (text: Bacteria 98.64, Archaea 1.36)",
   paste(rownames(dom_counts), round(100 * rowSums(dom_counts) / sum(dom_counts), 2)))
cat("\nArchaeal counts per sample (expected: > 0 only in marine CF samples):\n")
if ("d__Archaea" %in% rownames(dom_counts)) print(dom_counts["d__Archaea", SAMPLE_ORDER])

fam_glom <- tax_glom(ps, taxrank = "Family")
tax_table(fam_glom)[, "Family"] <- gsub("^f__", "", tax_table(fam_glom)[, "Family"])
fam_counts <- as(otu_table(fam_glom), "matrix")
rownames(fam_counts) <- as.vector(tax_table(fam_glom)[, "Family"])
fam_counts <- fam_counts[, SAMPLE_ORDER]
FAM <- t(sweep(fam_counts, 2, colSums(fam_counts), "/"))      

fam_pooled <- sort(100 * rowSums(fam_counts) / sum(depth$Depth), decreasing = TRUE)
ms("3.2 number of families (text: 123)", nrow(fam_counts))
ms("3.2 top-3 families, % of all counts (text: 22.05, 12.62, 9.73)",
   paste(names(fam_pooled)[1:3], round(fam_pooled[1:3], 2)))

# ---- A4. Figure 2A-B: families -------------------------------------------------
fam_long <- psmelt(fam_glom) %>% as_tibble() %>%
  mutate(Family = if_else(str_detect(Family, "^Incongruent|^Not_annotated"),
                          "Unclassified", Family))
top10_fam <- fam_long %>%
  filter(Family != "Unclassified") %>%
  group_by(Family) %>% summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>% slice_head(n = 10) %>% pull(Family)
fam_levels  <- c(top10_fam, "Others", "Unclassified")
fam_palette <- c(setNames(BASE_COLS[seq_along(top10_fam)], top10_fam),
                 Others = COL_OTHERS, Unclassified = COL_UNCLASS)

fam_plot_df <- fam_long %>%
  mutate(Family_plot = factor(if_else(Family %in% c(top10_fam, "Unclassified"),
                                      Family, "Others"), levels = fam_levels),
         Sample = factor(Sample, levels = SAMPLE_ORDER)) %>%
  group_by(Sample, Family_plot) %>% summarise(Abundance = sum(Abundance))

p_abs <- ggplot(fam_plot_df, aes(Sample, Abundance, fill = Family_plot)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = fam_palette, name = "Family", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scales::comma) +
  labs(x = LAB$stations, y = LAB$abs_tax) +
  theme_nature() + theme(legend.position = "none")

p_rel <- ggplot(fam_plot_df, aes(Sample, Abundance, fill = Family_plot)) +
  geom_col(position = "fill", width = 0.85) +
  scale_fill_manual(values = fam_palette, name = "Family", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = LAB$stations, y = LAB$rel) +
  theme_nature() + theme(legend.position = "right")

target_fams <- intersect(c("D2472", "Pelagibacteraceae", "Nanopelagicaceae",
                           "Methylophilaceae"), colnames(FAM))
cat("\nSelected families, % of family-assigned counts per sample:\n")
print(round(100 * FAM[, target_fams, drop = FALSE], 1))

# ---- A5. Table 1 and Figure S1: phyla ----------------------------------------
phy_glom <- tax_glom(ps, taxrank = "Phylum")
tax_table(phy_glom)[, "Phylum"] <- gsub("^p__", "", tax_table(phy_glom)[, "Phylum"])
phy_counts <- as(otu_table(phy_glom), "matrix")

table1 <- tibble(Phylum = as.vector(tax_table(phy_glom)[, "Phylum"]),
                 CF_reads = rowSums(phy_counts[, CF_LAB, drop = FALSE]),
                 VF_reads = rowSums(phy_counts[, VF_LAB, drop = FALSE])) %>%
  mutate(CF_pct = round(100 * CF_reads / sum(CF_reads), 2),
         VF_pct = round(100 * VF_reads / sum(VF_reads), 2),
         Total_pct = round(100 * (CF_reads + VF_reads) / sum(CF_reads + VF_reads), 2)) %>%
  arrange(desc(Total_pct))
print(table1, n = Inf)
save_tab(table1, "Table1_phylum_composition")

n_top_phyla <- 7
phy_long <- psmelt(phy_glom) %>% as_tibble() %>%
  mutate(Phylum = if_else(str_detect(Phylum, "^Incongruent|^Not_annotated"),
                          "Unclassified", Phylum))
top_phy <- phy_long %>% filter(Phylum != "Unclassified") %>%
  group_by(Phylum) %>% summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total)) %>% slice_head(n = n_top_phyla) %>% pull(Phylum)
phy_levels  <- c(top_phy, "Others", "Unclassified")
phy_palette <- c(setNames(BASE_COLS[seq_along(top_phy)], top_phy),
                 Others = COL_OTHERS, Unclassified = COL_UNCLASS)
phy_plot_df <- phy_long %>%
  mutate(Phylum_plot = factor(if_else(Phylum %in% c(top_phy, "Unclassified"),
                                      Phylum, "Others"), levels = phy_levels),
         Sample = factor(Sample, levels = SAMPLE_ORDER)) %>%
  group_by(Sample, Phylum_plot) %>% summarise(Abundance = sum(Abundance))

p_abs_phy <- ggplot(phy_plot_df, aes(Sample, Abundance, fill = Phylum_plot)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = phy_palette, name = "Phylum", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scales::comma) +
  labs(x = LAB$stations, y = LAB$abs_tax) +
  theme_nature() + theme(legend.position = "none")
p_rel_phy <- ggplot(phy_plot_df, aes(Sample, Abundance, fill = Phylum_plot)) +
  geom_col(position = "fill", width = 0.85) +
  scale_fill_manual(values = phy_palette, name = "Phylum", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = LAB$stations, y = LAB$rel) +
  theme_nature() + theme(legend.position = "right")
figS1 <- (p_abs_phy + p_rel_phy) + plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(figS1, "FigureS1_phylum_composition", 183, 80)

# ---- A6. Detection limit------------------------------------------------------
marker_fams <- intersect(c("D2472", "Marinisomataceae", "Thioglobaceae"), rownames(fam_counts))
detection <- map_dfr(seq_len(nrow(PAIRS)), function(i) {
  cf <- PAIRS$cf_lab[i]; vf <- PAIRS$vf_lab[i]
  n_vf <- sum(fam_counts[, vf]); N_cf <- sum(fam_counts[, cf])
  out <- tibble(Zone = PAIRS$zone[i], CF = cf, VF = vf, Taxon = marker_fams,
                CF_reads = fam_counts[marker_fams, cf], CF_depth = N_cf,
                VF_reads = fam_counts[marker_fams, vf], VF_depth = n_vf)
  if ("d__Archaea" %in% rownames(dom_counts))
    out <- bind_rows(out, tibble(Zone = PAIRS$zone[i], CF = cf, VF = vf, Taxon = "Archaea",
                                 CF_reads = dom_counts["d__Archaea", cf],
                                 CF_depth = sum(dom_counts[, cf]),
                                 VF_reads = dom_counts["d__Archaea", vf],
                                 VF_depth = sum(dom_counts[, vf])))
  out}) %>%
  mutate(CF_pct = round(100 * CF_reads / CF_depth, 2),
         expected_in_VF = round(VF_depth * CF_reads / CF_depth, 1),
         P_missed_by_chance = signif(dhyper(0, CF_reads, CF_depth - CF_reads,
                                            pmin(VF_depth, CF_depth)), 3))
cat("\nDetection limit (P_missed_by_chance = probability of zero reads in a\n",
    "random subsample of the CF sample at the VF depth):\n", sep = "")
print(detection, n = Inf, width = Inf)
save_tab(detection, "check_detection_limit")

for (i in seq_len(nrow(PAIRS))) {           
  r <- tryCatch({
    cf_rar <- rarefy_even_depth(prune_samples(PAIRS$cf_lab[i], fam_glom),
                                sample.size = sum(fam_counts[, PAIRS$vf_lab[i]]),
                                rngseed = CFG$seed, replace = FALSE, verbose = FALSE)
    fams <- as.vector(tax_table(cf_rar)[taxa_sums(cf_rar) > 0, "Family"])
    paste(marker_fams, ifelse(marker_fams %in% fams, "yes", "no"), collapse = " | ")
  }, error = function(e) paste("rarefy_even_depth failed:", conditionMessage(e)))
  cat(sprintf("  Station %s rarefied to its VF depth: %s\n", PAIRS$cf_lab[i], r))
}

# ---- A7. Alpha diversity (Figure 2C) --------------------------------------------
shannon <- estimate_richness(ps, measures = "Shannon") %>%
  rownames_to_column("Label") %>% mutate(Label = gsub("^X", "", Label)) %>%
  left_join(depth, by = "Label")
sp <- suppressWarnings(cor.test(shannon$Shannon, log10(shannon$Depth), method = "spearman"))
ms("2.4/3.2 Shannon vs depth, Spearman rho and P (text: 0.59, 0.029)",
   c(round(unname(sp$estimate), 2), signif(sp$p.value, 2)))

counts_mat <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) counts_mat <- t(counts_mat)
counts_mat <- counts_mat[, colSums(counts_mat) > 0, drop = FALSE]
if (any(counts_mat != round(counts_mat))) {
  warning("Non-integer mOTU counts were rounded for iNEXT.")
  counts_mat <- round(counts_mat)}

inext_in <- as.data.frame(counts_mat)
set.seed(CFG$seed)
cov_min <- min(iNEXT::iNEXT(inext_in, q = 1, datatype = "abundance")$DataInfo$SC)
set.seed(CFG$seed)
hill_q1 <- iNEXT::estimateD(inext_in, q = 1, datatype = "abundance",
                            base = "coverage", level = cov_min) %>%
  transmute(Label = gsub("^X", "", Assemblage), Hill_q1 = qD, LCL = qD.LCL, UCL = qD.UCL) %>%
  left_join(depth, by = "Label")
ms("2.4 common sample coverage used for Hill q1", round(cov_min, 4))
print(hill_q1 %>% dplyr::select(Label, Fraction, Zone, Hill_q1, LCL, UCL))
save_tab(hill_q1, "alpha_taxonomic_hill_q1")

cf_hill <- hill_q1 %>% filter(Fraction == "Cellular") %>% mutate(Zone = droplevels(Zone))
kw_hill <- kruskal.test(Hill_q1 ~ Zone, data = cf_hill)
kw_shan <- kruskal.test(Shannon ~ Zone, data = shannon %>% filter(Fraction == "Cellular"))
ms("3.2 Kruskal-Wallis Hill q1 ~ zone, CF: chi2, df, P (text: 5.95, 3, 0.114)",
   c(round(unname(kw_hill$statistic), 2), unname(kw_hill$parameter), signif(kw_hill$p.value, 3)))
ms("    (same test on raw Shannon, for comparison)",
   c(round(unname(kw_shan$statistic), 2), unname(kw_shan$parameter), signif(kw_shan$p.value, 3)))

if (CFG$run_new) {
  w <- wilcox.test(Hill_q1 ~ Group, data = cf_hill)
  ms("[NEW] Wilcoxon Hill q1, marine-influenced vs low-salinity CF: W, P",
     c(unname(w$statistic), signif(w$p.value, 3)))}

p_hill <- ggplot(hill_q1 %>% mutate(Label = factor(Label, levels = SAMPLE_ORDER)),
                 aes(Label, Hill_q1, colour = Fraction, shape = Fraction)) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0, linewidth = 0.4, alpha = 0.6) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = PAL$frac_tax, name = "Fraction") +
  scale_shape_manual(values = FRAC_SHAPES, name = "Fraction") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = LAB$stations, y = LAB$hill_tax) +
  theme_nature() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

# ---- A8. Beta diversity (Figure 2D), PERMANOVA and PERMDISP ------------------
ps_rel   <- transform_sample_counts(ps, function(x) x / sum(x))
bray_tax <- phyloseq::distance(ps_rel, method = "bray")
pcoa_tax <- ordinate(ps_rel, method = "PCoA", distance = bray_tax)
eig_tax  <- round(100 * pcoa_tax$values$Relative_eig[1:2], 1)
ms("3.2 PCoA variance, axes 1-2 (text: 67.8, 15.0)", eig_tax)

pcoa_tax_df <- as.data.frame(pcoa_tax$vectors[, 1:2]) %>%
  rownames_to_column("Label") %>% left_join(DESIGN_LAB, by = "Label")
p_pcoa_tax <- ggplot(pcoa_tax_df, aes(Axis.1, Axis.2, colour = Zone, shape = Fraction)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = PAL$zone_col, name = "Zone") +
  scale_shape_manual(values = FRAC_SHAPES, name = "Fraction") +
  labs(x = paste0("PCoA1 (", eig_tax[1], "%)"), y = paste0("PCoA2 (", eig_tax[2], "%)")) +
  theme_nature() + theme(legend.position = "right")

# PERMANOVA on the CF only (n = 10), two groups defined by measured salinity
bray_cf <- as.dist(as.matrix(bray_tax)[CF_LAB, CF_LAB])
meta_cf <- as.data.frame(DESIGN_LAB[match(CF_LAB, DESIGN_LAB$Label), ])
rownames(meta_cf) <- meta_cf$Label
set.seed(CFG$seed)
perm_tax <- adonis2(bray_cf ~ Group, data = meta_cf, permutations = 999)
set.seed(CFG$seed)
disp_tax <- permutest(betadisper(bray_cf, meta_cf$Group), permutations = 999)
print(perm_tax)
ms("3.2 PERMANOVA CF, two groups: F, R2, P (text: 42.9, 0.84, 0.005)",
   c(round(perm_tax$F[1], 1), round(perm_tax$R2[1], 2), perm_tax$`Pr(>F)`[1]))
ms("3.2 PERMDISP P (text: 0.67)", disp_tax$tab$`Pr(>F)`[1])
if (CFG$run_new) {
  ex <- exact_p_two_groups(bray_cf, meta_cf$Group)
  ms("[NEW] exact P by enumeration: F, allocations, P", round(ex, 4))
  # If P_exact = 1/210, write "P = 0.005, the minimum attainable with 4 vs 6 samples".
}

# Robustness: CLR + Aitchison distance 
if (has_pkg("microbiome")) {
  ps_clr   <- microbiome::transform(ps, "clr")
  pcoa_ait <- ordinate(ps_clr, method = "PCoA",
                       distance = phyloseq::distance(ps_clr, method = "euclidean"))
  set.seed(CFG$seed)
  pt <- protest(pcoa_tax$vectors[, 1:2], pcoa_ait$vectors[, 1:2], permutations = 999)
  ms("robustness: Procrustes correlation Bray vs Aitchison, P",
     c(round(sqrt(1 - pt$ss), 2), pt$signif))
}

# Paired CF-VF Bray-Curtis dissimilarity at family level (descriptive, n = 4).
paired_bc <- PAIRS %>%
  mutate(BC_family = map2_dbl(cf_lab, vf_lab,
                              ~ as.numeric(vegdist(FAM[c(.x, .y), ], method = "bray")))) %>%
  dplyr::select(zone, cf_lab, vf_lab, BC_family)
cat("\nPaired CF-VF Bray-Curtis dissimilarity (family level):\n"); print(paired_bc)
save_tab(paired_bc, "check_paired_BC_taxonomy")

# ---- A9. Sensitivity to typical reagent contaminants ----------------------
if (CFG$run_new) {
  suspects <- c("Cutibacterium", "Stenotrophomonas", "Pseudomonas", "Brevundimonas",
                "Sphingomonas", "Acinetobacter", "Ralstonia", "Bradyrhizobium",
                "Methylobacterium", "Staphylococcus", "Corynebacterium", "Acidovorax")
  genus <- str_remove(gsub("^g__", "", as.vector(tax_table(ps)[, "Genus"])), "_[A-Z]+$")
  is_suspect <- !is.na(genus) & genus %in% suspects
  cnt <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) cnt <- t(cnt)
  suspect_share <- tibble(Label = colnames(cnt),
                          suspect_pct = round(100 * colSums(cnt[is_suspect, , drop = FALSE]) /
                                                colSums(cnt), 2)) %>%
    left_join(DESIGN_LAB, by = "Label") %>% dplyr::select(Label, Fraction, Zone, suspect_pct)
  cat("\n[NEW] Share of counts in typical contaminant genera (%):\n"); print(suspect_share, n = Inf)
  if (any(is_suspect) && sum(!is_suspect) > 1) {
    ps_clean <- prune_taxa(taxa_names(ps)[!is_suspect], ps)
    rel_clean <- transform_sample_counts(ps_clean, function(x) x / sum(x))
    bc_clean  <- as.matrix(phyloseq::distance(rel_clean, method = "bray"))
    bc_full   <- as.matrix(bray_tax)
    sens <- PAIRS %>% transmute(zone, pair = paste(cf_lab, vf_lab, sep = " - "),
                                BC_all_taxa = map2_dbl(cf_lab, vf_lab, ~ bc_full[.x, .y]),
                                BC_without_suspects = map2_dbl(cf_lab, vf_lab, ~ bc_clean[.x, .y]))
    cat("\n[NEW] Paired CF-VF Bray-Curtis (mOTU level) with and without them:\n"); print(sens)
    save_tab(sens, "check_contaminant_sensitivity")}}

# ---- A10. Figure 2 ----------------------------------------------------------------
fig2 <- (p_abs + p_rel) / (p_hill + p_pcoa_tax) +
  plot_layout(heights = c(3, 2)) + plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig2, "Figure2", 183, 150)


# ==============================================================================
# PART B. FUNCTION (KEGG orthologues) -- Figure 3, Section 3.3
# ==============================================================================

# ---- B1. Import SqueezeMeta tables ---------------------------------------------
fix_z3_2 <- function(df) rename_with(df, ~ str_replace(.x, "Z3_2", "Z4_2"))

kegg <- read_tsv(CFG$in_kegg_tpm, show_col_types = FALSE) %>%
  fix_z3_2() %>%
  mutate(KEGG_ID = str_trim(str_remove(KEGG_ID, "^ko:"))) %>%
  filter(!is.na(KEGG_ID), KEGG_ID != "",
         !str_detect(KEGG_ID, ";"))          # annotations with several KOs are dropped
raw_cols <- names(kegg)[str_starts(names(kegg), "Raw_")]
tpm_cols <- names(kegg)[str_starts(names(kegg), "TPM_")]

orf_tax <- read_tsv(CFG$in_kegg_tax, show_col_types = FALSE)
if ("Estacion" %in% names(orf_tax)) orf_tax <- dplyr::rename(orf_tax, Station = Estacion)
orf_tax <- orf_tax %>%
  mutate(Station = str_replace(Station, "^Z3_2$", "Z4_2"),
         KEGG_ID = str_trim(str_remove(KEGG_ID, "^ko:")),
         across(any_of(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")),
                ~ if_else(is.na(.x) | .x == "" | .x == "No definido", "Unassigned", .x))) %>%
  filter(Station %in% DESIGN$sqm_id)

env <- read_excel(CFG$in_env) %>% rename_with(trimws)

# ---- B2. KO x station tables -------------------------------------------------
long_from <- function(cols, prefix) {
  kegg %>% dplyr::select(KEGG_ID, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "Station", values_to = "value") %>%
    mutate(Station = str_remove(Station, prefix)) %>%
    filter(Station %in% DESIGN$sqm_id, !is.na(value))
}
raw_rows <- long_from(raw_cols, "^Raw_")
tpm_rows <- long_from(tpm_cols, "^TPM_")

dup <- raw_rows %>% filter(value > 0) %>% dplyr::count(KEGG_ID, Station) %>% filter(n > 1)
if (nrow(dup)) warning(nrow(dup), " KO x station combinations have counts in more than one ",
                       "row of kegg_tpm.tsv; they are summed here.")

tpm_first <- tpm_rows %>% distinct(KEGG_ID, Station, .keep_all = TRUE)
tpm_sum   <- tpm_rows %>% group_by(KEGG_ID, Station) %>% summarise(value = sum(value))
n_diff <- inner_join(tpm_first, tpm_sum, by = c("KEGG_ID", "Station"),
                     suffix = c("_first", "_sum")) %>%
  filter(abs(value_first - value_sum) > 1e-9) %>% nrow()
cat("\n[CHECK] KO x station TPM values that differ between 'first row' and 'sum':", n_diff, "\n")

raw_long <- raw_rows %>% group_by(KEGG_ID, Station) %>% summarise(value = sum(value)) %>%
  filter(value > 0)
tpm_long <- tpm_sum %>% filter(value > 0)

comm_raw <- raw_long %>%
  mutate(Station = factor(Station, levels = DESIGN$sqm_id)) %>%
  group_by(Station, KEGG_ID) %>% summarise(Raw = sum(value)) %>%
  pivot_wider(names_from = KEGG_ID, values_from = Raw, values_fill = 0) %>%
  mutate(Station = as.character(Station)) %>%
  column_to_rownames("Station") %>% as.matrix() %>% round()
stopifnot(setequal(rownames(comm_raw), DESIGN$sqm_id))
counts_ko <- t(comm_raw)[, DESIGN$sqm_id]                    # KO x station (integers)
tpm_ko <- tpm_long %>%
  pivot_wider(names_from = Station, values_from = value, values_fill = 0) %>%
  column_to_rownames("KEGG_ID") %>% as.matrix()
tpm_ko <- tpm_ko[, DESIGN$sqm_id]

ko_annot <- kegg %>% group_by(KEGG_ID) %>%
  summarise(across(any_of(c("gene_name", "KEGGFUN", "KEGGPATH")), ~ dplyr::first(na.omit(.x))))

ms("3.3 KOs detected (Fig. 3 legend: 8,829)", nrow(counts_ko))
ko_depth <- tibble(Station = colnames(counts_ko), Depth = colSums(counts_ko)) %>%
  left_join(DESIGN_SQM, by = "Station")
ms("KO-assigned reads, median CF / VF",
   c(median(ko_depth$Depth[ko_depth$Fraction == "Cellular"]),
     median(ko_depth$Depth[ko_depth$Fraction == "Vesicular"])))

# ---- B3. KEGG level-2 pathway categories (Figure 3A-B) ---------------------------
PATHWAY_L1 <- c("Metabolism", "Genetic Information Processing",
                "Environmental Information Processing", "Cellular Processes",
                "Organismal Systems", "Human Diseases")
THRESH <- 1.0      # categories below this mean relative abundance (%) -> "Others"

path_tbl <- ko_annot %>%
  dplyr::select(KEGG_ID, KEGGPATH) %>%
  filter(!is.na(KEGGPATH), str_trim(KEGGPATH) != "") %>%
  separate_rows(KEGGPATH, sep = "\\s*\\|\\s*") %>%
  mutate(KEGGPATH = str_trim(KEGGPATH)) %>% filter(KEGGPATH != "") %>%
  separate(KEGGPATH, into = c("L1", "L2", "L3"), sep = "; ", fill = "right", extra = "merge") %>%
  mutate(across(c(L1, L2, L3), str_trim)) %>%
  filter(L1 %in% PATHWAY_L1, !is.na(L2), L2 != "") %>%
  distinct(KEGG_ID, L2) %>%
  group_by(KEGG_ID) %>% summarise(L2 = dplyr::first(L2))     # first-listed category
ms("3.3 KOs with a pathway category: n, % (legend: 4,973, 56%)",
   c(nrow(path_tbl), round(100 * nrow(path_tbl) / nrow(counts_ko))))

agg_cat <- function(d) {
  d %>% inner_join(path_tbl, by = "KEGG_ID") %>%
    group_by(Station, category = L2) %>% summarise(val = sum(value)) %>% ungroup()}

raw_cat <- agg_cat(raw_long); tpm_cat <- agg_cat(tpm_long)
keep_cats <- tpm_cat %>% group_by(Station) %>% mutate(p = 100 * val / sum(val)) %>%
  group_by(category) %>% summarise(mean_pct = mean(p)) %>%
  filter(mean_pct >= THRESH) %>% pull(category)

collapse_cat <- function(d) {
  d %>% mutate(category = if_else(category %in% keep_cats, category, "Others")) %>%
    group_by(Station, category) %>% summarise(val = sum(val)) %>% ungroup()}

raw_cat <- collapse_cat(raw_cat) %>% left_join(DESIGN_SQM, by = "Station")
tpm_cat <- collapse_cat(tpm_cat) %>% group_by(Station) %>%
  mutate(rel_pct = 100 * val / sum(val)) %>% ungroup() %>% left_join(DESIGN_SQM, by = "Station")

cat_by_fraction <- tpm_cat %>% group_by(Fraction, category) %>%
  summarise(mean_pct = round(mean(rel_pct), 1)) %>%
  pivot_wider(names_from = Fraction, values_from = mean_pct, values_fill = 0) %>%
  arrange(desc(Cellular))

print(cat_by_fraction, n = Inf)
save_tab(cat_by_fraction, "TableS_KEGG_categories_by_fraction")

cat_order <- c(tpm_cat %>% filter(Fraction == "Cellular", category != "Others") %>%
                 group_by(category) %>% summarise(t = sum(val)) %>% arrange(t) %>% pull(category),
               "Others")
named_cats <- setdiff(cat_order, "Others")

kegg_palette <- c(setNames(rev(BASE_COLS[seq_along(named_cats)]), named_cats),
                  Others = COL_UNCLASS)
lab_sqm <- setNames(DESIGN$label, DESIGN$sqm_id)
mk_bar <- function(d, yvar, ylab, scale = 1) {
  d %>% mutate(category = factor(category, levels = cat_order),
               Station = factor(Station, levels = DESIGN$sqm_id),
               y = .data[[yvar]] / scale) %>%
    ggplot(aes(Station, y, fill = category)) +
    geom_col(width = 0.9, colour = "black", linewidth = 0.15) +
    scale_fill_manual(values = kegg_palette, drop = FALSE, name = LAB$kegg_cat) +
    scale_x_discrete(labels = lab_sqm) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = LAB$stations, y = ylab) + theme_nature()
}

p_ko_abs <- mk_bar(raw_cat, "val", LAB$abs_ko, scale = 1e6)
p_ko_rel <- mk_bar(tpm_cat, "rel_pct", LAB$rel)

# ---- B4. KO alpha diversity (Figure 3C-D) and supplementary table ------------
alpha_ko <- tibble(Station = colnames(counts_ko),
                   Depth = colSums(counts_ko),
                   q0 = colSums(counts_ko > 0),
                   q1 = exp(vegan::diversity(t(counts_ko), index = "shannon")),
                   Coverage = apply(counts_ko, 2, good_turing)) %>%
  left_join(DESIGN_SQM, by = "Station")
print(alpha_ko %>% dplyr::select(Label, Fraction, Zone, Depth, q0, q1, Coverage) %>%
        mutate(q1 = round(q1, 1), Coverage = signif(Coverage, 8)), n = Inf)
save_tab(alpha_ko, "TableS_KO_alpha_diversity_coverage")

w0 <- wilcox.test(q0 ~ Fraction, data = alpha_ko)
w1 <- wilcox.test(q1 ~ Fraction, data = alpha_ko)
 
mk_alpha <- function(yvar, ylab) {
  alpha_ko %>% mutate(Station = factor(Station, levels = DESIGN$sqm_id)) %>%
    ggplot(aes(Station, .data[[yvar]], colour = Zone, shape = Fraction)) +
    geom_point(size = 3, alpha = 0.9) +
    scale_colour_manual(values = PAL$zone_col, name = "Zone", limits = ZONE_LEVELS) +
    scale_shape_manual(values = FRAC_SHAPES, name = "Fraction", limits = FRAC_LEVELS) +
    scale_x_discrete(labels = lab_sqm) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = LAB$stations, y = ylab) + theme_nature()}

p_q0 <- mk_alpha("q0", LAB$ko_rich)
p_q1 <- mk_alpha("q1", LAB$ko_hill)

# ---- B5. KO beta diversity (Figure 3E), PERMANOVA and PERMDISP ---------------------
comm_tpm <- t(tpm_ko)
bc_fun <- vegdist(comm_tpm, method = "bray")
md <- as.data.frame(DESIGN_SQM); rownames(md) <- md$Station
md <- md[rownames(comm_tpm), ]

# (i) fraction, all 14 samples
set.seed(CFG$seed); pm_frac <- adonis2(bc_fun ~ Fraction, data = md, permutations = 9999)
set.seed(CFG$seed); disp_frac <- permutest(betadisper(bc_fun, md$Fraction), permutations = 999)
# (ii) zone within the CF, without the unreplicated Coastal Waters station (n = 9)
keep9 <- rownames(md)[md$Fraction == "Cellular" & md$Zone != "Coastal Waters"]
bc_cf9 <- as.dist(as.matrix(bc_fun)[keep9, keep9]); md9 <- droplevels(md[keep9, ])
set.seed(CFG$seed); pm_zone   <- adonis2(bc_cf9 ~ Zone, data = md9, permutations = 9999)
set.seed(CFG$seed); disp_zone <- permutest(betadisper(bc_cf9, md9$Zone), permutations = 999)

pcoa_fun <- cmdscale(bc_fun, k = 2, eig = TRUE)
eig_fun  <- round(100 * pcoa_fun$eig[1:2] / sum(pcoa_fun$eig[pcoa_fun$eig > 0]), 1)
ms("3.3 PCoA variance, axes 1-2 (text: 73.8, 11.8)", eig_fun)
 
pcoa_fun_df <- tibble(Station = rownames(pcoa_fun$points),
                      Axis1 = pcoa_fun$points[, 1], Axis2 = pcoa_fun$points[, 2]) %>%
  left_join(DESIGN_SQM, by = "Station")
p_pcoa_fun <- ggplot(pcoa_fun_df, aes(Axis1, Axis2, colour = Zone, shape = Fraction)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_colour_manual(values = PAL$zone_col, name = "Zone", limits = ZONE_LEVELS) +
  scale_shape_manual(values = FRAC_SHAPES, name = "Fraction", limits = FRAC_LEVELS) +
  labs(x = paste0("PCoA1 (", eig_fun[1], "%)"), y = paste0("PCoA2 (", eig_fun[2], "%)")) +
  theme_nature()

# ---- B6. Figure 3 --------------------------------------------------------------------
row1 <- (p_ko_abs + p_ko_rel) + plot_layout(guides = "collect") & theme(legend.position = "right")
row2 <- (p_q0 + p_q1 + p_pcoa_fun) + plot_layout(guides = "collect") & theme(legend.position = "right")
fig3 <- row1 / row2 + plot_layout(heights = c(3, 2)) + plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig3, "Figure3", 260, 170)


# ==============================================================================
# PART C. CF-VF PARTITION OF KOs -- Figure 4, Section 3.4
# ==============================================================================

# ---- C1. KEGG pathway sets ---------------------------------------------------------
if (all(file.exists(c(CFG$in_ko_path, CFG$in_path_names)))) {
  ko_path_tbl <- readRDS(CFG$in_ko_path)
  pathway_names_tbl <- readRDS(CFG$in_path_names)
} else if (has_pkg("KEGGREST")) {
  message("KEGG tables not found: downloading with KEGGREST (", Sys.Date(), ").")
  lk <- KEGGREST::keggLink("pathway", "ko")
  pn <- KEGGREST::keggList("pathway")
  ko_path_tbl <- tibble(KEGG_ID = str_remove(names(lk), "^ko:"),
                        Pathway_ID = str_remove(unname(lk), "^path:")) %>%
    filter(str_starts(Pathway_ID, "map"))
  pathway_names_tbl <- tibble(Pathway_ID = str_remove(names(pn), "^path:"), Pathway = unname(pn))
  saveRDS(ko_path_tbl, CFG$in_ko_path); saveRDS(pathway_names_tbl, CFG$in_path_names)
} else stop("Provide ko_path_tbl.rds and pathway_names_tbl.rds, or install KEGGREST.")

ko_names <- if (file.exists(CFG$in_ko_info)) {
  readRDS(CFG$in_ko_info) %>%
    transmute(KEGG_ID, name = coalesce(Symbol, Protein, KEGG_ID)) %>%
    distinct(KEGG_ID, .keep_all = TRUE)
} else if ("gene_name" %in% names(ko_annot)) {
  ko_annot %>% transmute(KEGG_ID, name = coalesce(gene_name, KEGG_ID))
} else tibble(KEGG_ID = rownames(counts_ko), name = rownames(counts_ko))

# Plant, animal and other eukaryote-specific maps receive hits from promiscuous
# enzyme families in a prokaryotic metagenome and are excluded from pathway
# statistics (list them in a supplementary table).
MAPS_EXCLUDED <- c("map00100", "map00140", "map00510", "map00513", "map00512",
                   "map00514", "map00515", "map00532", "map00534", "map00533",
                   "map00531", "map00601", "map00603", "map00604", "map00590",
                   "map00591", "map00592", "map00565", "map00941", "map00940")

TERM2GENE <- ko_path_tbl %>%
  distinct(Pathway_ID, KEGG_ID) %>%
  filter(str_starts(Pathway_ID, "map00") | str_starts(Pathway_ID, "map02"),   
         !Pathway_ID %in% (MAPS_EXCLUDED),              
         KEGG_ID %in% rownames(counts_ko)) %>%
  as.data.frame()
TERM2NAME <- pathway_names_tbl %>%
  distinct(Pathway_ID, Pathway) %>%
  filter(Pathway_ID %in% TERM2GENE$Pathway_ID) %>%
  right_join(tibble(Pathway_ID = unique(TERM2GENE$Pathway_ID)), by = "Pathway_ID") %>%
  mutate(Pathway = coalesce(Pathway, Pathway_ID)) %>%
  as.data.frame()
ms("2.6 pathway maps used / KOs with a map", c(n_distinct(TERM2GENE$Pathway_ID),
                                               n_distinct(TERM2GENE$KEGG_ID)))

# ---- C2. Rarefaction and KO partition (Figure 4A and 4C) ---------------------
comm_pairs <- comm_raw[rownames(comm_raw) %in% c(PAIRS$cf, PAIRS$vf), ]
rar_depth  <- min(rowSums(comm_pairs))
ms("3.4 rarefaction depth (text: 1,815,572)", rar_depth, digits = 9)
set.seed(CFG$seed)
comm_rar <- rrarefy(comm_pairs, sample = rar_depth)      # one draw, as in the source
pa <- comm_rar > 0
kos_in <- function(st) colnames(pa)[pa[st, ]]

zone_part <- map_dfr(seq_len(nrow(PAIRS)), function(i) {
  cf <- kos_in(PAIRS$cf[i]); vf <- kos_in(PAIRS$vf[i])
  tibble(zone = PAIRS$zone[i], richness_CF = length(cf), richness_VF = length(vf),
         total = length(union(cf, vf)), shared = length(intersect(cf, vf)),
         only_CF = length(setdiff(cf, vf)), only_VF = length(setdiff(vf, cf)),
         VF_pct_of_CF = round(100 * length(vf) / length(cf), 1))
})

print(zone_part)
save_tab(zone_part, "TableS_KO_partition_by_zone")

kos_cf  <- colnames(pa)[colSums(pa[PAIRS$cf, , drop = FALSE]) > 0]
kos_vf  <- colnames(pa)[colSums(pa[PAIRS$vf, , drop = FALSE]) > 0]
excl_cf <- setdiff(kos_cf, kos_vf); excl_vf <- setdiff(kos_vf, kos_cf)
shared  <- intersect(kos_cf, kos_vf)

p_line <- ggplot(zone_part, aes(zone, total, group = 1)) +
  geom_line(colour = "grey50", linewidth = 0.6) +
  geom_point(size = 2.5, colour = "#3C8DAD") +
  geom_text(aes(label = total), vjust = -1, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  labs(x = NULL, y = LAB$ko_rich) + theme_tfm()

if (has_pkg("eulerr")) {
  euler_panels <- pmap(list(as.character(zone_part$zone), zone_part$only_CF,
                            zone_part$only_VF, zone_part$shared),
                       function(z, a, b, ab) {
                         fit <- eulerr::euler(c("Cellular" = a, "Vesicular" = b,
                                                "Cellular&Vesicular" = ab))
                         wrap_elements(full = plot(
                           fit, fills = list(fill = unname(PAL$frac_col), alpha = 0.6),
                           edges = list(col = "grey30"), labels = list(fontsize = 6),
                           quantities = list(fontsize = 6), main = list(label = z, fontsize = 7)))
                       })
  p_partition <- wrap_plots(euler_panels, nrow = 1) / p_line + plot_layout(heights = c(1, 1.1))
} else {
  p_partition <- p_line
}

pooled_df <- tibble(group = factor(c("Cellular only", "Vesicular only", "Shared"),
                                   levels = c("Shared", "Vesicular only", "Cellular only")),
                    n = c(length(excl_cf), length(excl_vf), length(shared)))
p_pooled <- ggplot(pooled_df, aes(n, group, fill = group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.25, size = 2.6, fontface = "bold") +
  scale_fill_manual(values = c("Cellular only" = unname(PAL$frac_col["Cellular"]),
                               "Vesicular only" = unname(PAL$frac_col["Vesicular"]),
                               "Shared" = "#9E9E9E"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = LAB$ko_n, y = NULL) + theme_tfm()

# ---- C3. ORFs with a KO (Figure 4B) --------------------------------------------------
orf_station <- orf_tax %>% dplyr::count(Station, name = "ORFs") %>%
  left_join(DESIGN_SQM, by = "Station") %>%
  left_join(ko_depth %>% dplyr::select(Station, Depth), by = "Station") %>%
  mutate(ORFs_per_M_reads = round(ORFs / (Depth / 1e6)))
print(orf_station %>% dplyr::select(Label, Fraction, Zone, ORFs, Depth, ORFs_per_M_reads), n = Inf)
ms("3.4 ORFs with KO, Inlet CF stations (text: 155,595-182,123)",
   range(orf_station$ORFs[orf_station$Fraction == "Cellular" & orf_station$Zone == "Inlet"]), digits = 7)
orf_pairs <- PAIRS %>%
  mutate(ORFs_CF = orf_station$ORFs[match(cf, orf_station$Station)],
         ORFs_VF = orf_station$ORFs[match(vf, orf_station$Station)],
         ratio = round(ORFs_CF / ORFs_VF, 1)) %>%
  dplyr::select(zone, ORFs_CF, ORFs_VF, ratio)

orf_zone <- orf_station %>% group_by(Zone, Fraction) %>% summarise(ORFs = mean(ORFs)) %>% ungroup()
p_orf <- ggplot(orf_zone, aes(Zone, ORFs, colour = Fraction, group = Fraction)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.5) +
  scale_colour_manual(values = PAL$frac_col, name = "Fraction") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = LAB$orf_n) + theme_tfm()

# ---- C4. Paired VF/CF ratios per pathway (Figure 4D) ---------------------------------
ko_prop <- tpm_sum %>%
  filter(Station %in% c(PAIRS$cf, PAIRS$vf), KEGG_ID %in% TERM2GENE$KEGG_ID) %>%
  left_join(DESIGN_SQM, by = "Station") %>%
  group_by(Station) %>% mutate(prop = value / sum(value)) %>% ungroup()

path_prop <- ko_prop %>%
  inner_join(TERM2GENE, by = "KEGG_ID", relationship = "many-to-many") %>%
  group_by(StationNum, Fraction, Pathway_ID) %>% summarise(path_prop = sum(prop)) %>% ungroup()
eps_path <- min(path_prop$path_prop[path_prop$path_prop > 0]) / 2

state_of <- function(cf, vf) case_when(cf > 0 & vf > 0 ~ "shared", cf == 0 & vf > 0 ~ "VF only",
                                       cf > 0 & vf == 0 ~ "CF only", TRUE ~ "absent")
ratio_pair <- path_prop %>%
  pivot_wider(names_from = Fraction, values_from = path_prop, values_fill = 0) %>%
  mutate(state = state_of(Cellular, Vesicular),
         log2_ratio = log2((Vesicular + eps_path) / (Cellular + eps_path)))
 
ratio_summary <- ratio_pair %>%
  group_by(Pathway_ID) %>%
  summarise(n_shared = sum(state == "shared"), n_absent = sum(state == "absent"),
            log2_ratio_mean = mean(log2_ratio[state == "shared"])) %>%
  filter(n_shared >= CFG$min_pairs) %>%
  left_join(TERM2NAME, by = "Pathway_ID")

ratio_plot_df <- bind_rows(slice_max(ratio_summary, log2_ratio_mean, n = 10),
                           slice_min(ratio_summary, log2_ratio_mean, n = 10)) %>%
  distinct(Pathway_ID, .keep_all = TRUE) %>%
  mutate(Higher_in = if_else(log2_ratio_mean > 0, "Vesicular", "Cellular"),
         Pathway = fct_reorder(Pathway, log2_ratio_mean))

print(ratio_plot_df %>% arrange(desc(log2_ratio_mean)) %>%
        transmute(Pathway_ID, Pathway = str_trunc(as.character(Pathway), 48),
                  n_shared, n_absent, log2 = round(log2_ratio_mean, 2)), n = Inf)
save_tab(ratio_summary %>% arrange(desc(log2_ratio_mean)), "TableS_pathway_log2_VF_CF")

p_ratio <- ggplot(ratio_plot_df, aes(log2_ratio_mean, Pathway, colour = Higher_in)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = log2_ratio_mean, yend = Pathway), linewidth = 0.6) +
  geom_point(aes(size = n_shared)) +
  scale_colour_manual(values = PAL$frac_col, name = NULL,
                      labels = c(Cellular = "Higher in the cellular fraction",
                                 Vesicular = "Higher in the vesicular fraction")) +
  scale_size_continuous(name = "Paired stations with the\npathway in both fractions",
                        range = c(1.2, 3.5), breaks = 1:4, limits = c(1, 4)) +
  labs(x = LAB$ratio, y = NULL) +
  theme_tfm() + theme(legend.position = "top", legend.box = "vertical")

GENERIC_KOS <- c(K00973 = "rmlA", K01710 = "rmlB", K01790 = "rmlC", K00067 = "rmlD",
                 K01428 = "ureC", K01429 = "ureB", K01430 = "ureA", K14048 = "ureAB",
                 K00626 = "atoB (thiolase)", K00632 = "fadA (thiolase)",
                 K01692 = "enoyl-CoA hydratase", K00558 = "dcm (DNA methyltransferase)")
ko_wide <- ko_prop %>% dplyr::select(StationNum, Fraction, KEGG_ID, prop) %>%
  pivot_wider(names_from = Fraction, values_from = prop, values_fill = 0)
eps_ko <- min(ko_prop$prop[ko_prop$prop > 0]) / 2
leaders <- map_dfr(as.character(ratio_plot_df$Pathway_ID), function(pid) {
  st <- ratio_pair$StationNum[ratio_pair$Pathway_ID == pid & ratio_pair$state == "shared"]
  kos <- TERM2GENE$KEGG_ID[TERM2GENE$Pathway_ID == pid]
  d <- ko_wide %>% filter(KEGG_ID %in% kos, StationNum %in% st) %>%
    group_by(KEGG_ID) %>% summarise(diff = sum(Vesicular - Cellular) / length(st))
  if (!nrow(d)) return(NULL)
  sgn <- if (sum(d$diff) >= 0) 1 else -1
  contrib <- pmax(sgn * d$diff, 0)
  d %>% mutate(Pathway_ID = pid, higher_in = if_else(sgn > 0, "VF", "CF"),
               contribution_pct = round(if (sum(contrib) > 0) 100 * contrib / sum(contrib) else 0, 1),
               generic = unname(GENERIC_KOS[KEGG_ID])) %>%
    arrange(desc(contribution_pct)) %>% slice_head(n = 5)
}) %>%
  left_join(ko_names, by = "KEGG_ID") %>% left_join(TERM2NAME, by = "Pathway_ID") %>%
  dplyr::select(Pathway, higher_in, KEGG_ID, name, generic, contribution_pct)
cat("\nLeading KOs of each pathway in Fig. 4D (top 5 by contribution):\n")
print(leaders %>% mutate(Pathway = str_trunc(Pathway, 34), name = str_trunc(name, 28)), n = Inf)
save_tab(leaders, "TableS_leading_KOs_Fig4D")

# KO-level ratios
ratio_ko <- ko_wide %>%
  mutate(state = state_of(Cellular, Vesicular),
         log2_ratio = log2((Vesicular + eps_ko) / (Cellular + eps_ko))) %>%
  group_by(KEGG_ID) %>%
  summarise(n_shared = sum(state == "shared"),
            log2_ratio_mean = mean(log2_ratio[state == "shared"])) %>%
  filter(n_shared >= CFG$min_pairs) %>%
  left_join(ko_names, by = "KEGG_ID") %>% arrange(desc(log2_ratio_mean))
save_tab(ratio_ko, "TableS_KO_log2_VF_CF")

# ---- C5. Over-representation analysis of fraction-exclusive KOs --------------
if (has_pkg("clusterProfiler")) {
  universe <- intersect(union(kos_cf, kos_vf), TERM2GENE$KEGG_ID)
  run_ora <- function(genes) {
    e <- clusterProfiler::enricher(gene = intersect(genes, universe), universe = universe,
                                   TERM2GENE = TERM2GENE[, c("Pathway_ID", "KEGG_ID")],
                                   TERM2NAME = TERM2NAME[, c("Pathway_ID", "Pathway")],
                                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                   minGSSize = 5, maxGSSize = 500)
    if (is.null(e)) return(tibble())
    as_tibble(as.data.frame(e)) %>% arrange(p.adjust)
  }
  ora_cf <- run_ora(excl_cf); ora_vf <- run_ora(excl_vf)
  cat("\n[MS] 3.4 ORA, CF-exclusive KOs, BH-adjusted P < 0.05 (text: toluene, benzoate,\n",
      "aminobenzoate, ethylbenzene, polycyclic aromatic hydrocarbon degradation):\n", sep = "")
  if (nrow(ora_cf)) print(ora_cf %>% filter(p.adjust < 0.05) %>%
                            dplyr::select(ID, Description, GeneRatio, BgRatio, p.adjust), n = Inf)
  if (nrow(ora_vf)) ms("3.4 ORA, VF-exclusive KOs: smallest adjusted P (text: 0.15)",
                       signif(min(ora_vf$p.adjust), 2))
  save_tab(ora_cf, "TableS_ORA_CF_exclusive"); save_tab(ora_vf, "TableS_ORA_VF_exclusive")
} else {
  ora_cf <- tibble()
  message("clusterProfiler not installed: ORA skipped.")}

# ---- C6. Exclusive KOs -------------------------------------------------------
if ("Phylum" %in% names(orf_tax)) {
  phylum_share <- function(kos) {
    orf_tax %>% filter(Station %in% PAIRS$cf, KEGG_ID %in% kos) %>%
      dplyr::count(Phylum, sort = TRUE) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
      slice_head(n = 6)}
 
    if (nrow(ora_cf) && any(ora_cf$p.adjust < 0.05)) {
    sig_kos <- intersect(excl_cf, TERM2GENE$KEGG_ID[TERM2GENE$Pathway_ID %in%
                                                      ora_cf$ID[ora_cf$p.adjust < 0.05]])
    print(phylum_share(sig_kos))}}

tetx <- ko_names %>% filter(KEGG_ID %in% excl_vf, str_detect(name, regex("tetX", ignore_case = TRUE)))

if (CFG$run_new) {
  if ("Kingdom" %in% names(orf_tax)) {
    kingdom_share <- orf_tax %>% dplyr::count(Station, Kingdom) %>%
      group_by(Station) %>% mutate(pct = round(100 * n / sum(n), 2)) %>% ungroup() %>%
      left_join(DESIGN_SQM %>% dplyr::select(Station, Label), by = "Station") %>%
      dplyr::select(Label, Kingdom, pct) %>%
      pivot_wider(names_from = Kingdom, values_from = pct, values_fill = 0)
    save_tab(kingdom_share, "check_kingdom_of_ORFs")
  }
  # (b) Reagent contaminants: which genera carry the VF-exclusive KOs
  if ("Genus" %in% names(orf_tax)) {
    vf_only_genus <- orf_tax %>% filter(Station %in% PAIRS$vf, KEGG_ID %in% excl_vf) %>%
      left_join(DESIGN_SQM %>% dplyr::select(Station, Label), by = "Station") %>%
      dplyr::count(Label, Genus) %>% group_by(Label) %>%
      mutate(pct = round(100 * n / sum(n), 1)) %>% slice_max(n, n = 8, with_ties = FALSE) %>% ungroup()
    cat("\n[NEW] Genus of ORFs carrying VF-exclusive KOs (top 8 per sample):\n")
    print(vf_only_genus, n = Inf)
    save_tab(vf_only_genus, "check_genus_of_VF_exclusive_KOs")
  }
}

# ---- C7. Figure 4 ------------------------------------------------------------
fig4 <- (wrap_elements(full = p_partition) / p_orf / p_pooled + plot_layout(heights = c(2, 1.2, 1)) |
           p_ratio) +
  plot_layout(widths = c(1, 1.1)) + plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig4, "Figure4", 260, 170)

# ==============================================================================
# PART D. NITROGEN CYCLE (NCycDB) -- Figure 5, Section 3.5
# ==============================================================================

tidy_process <- function(x) str_to_sentence(str_squish(as.character(x)))

# ---- D1. Import --------------------------------------------------------------
ncyc_raw <- read_excel(CFG$in_ncyc, sheet = CFG$in_ncyc_sheet, skip = 1, col_names = TRUE)
names(ncyc_raw)[1:2] <- c("Pathway", "Symbol")
ncyc <- ncyc_raw %>%
  dplyr::select(Pathway, Symbol, ends_with("_TPM")) %>%
  pivot_longer(-c(Pathway, Symbol), names_to = "Station", values_to = "TPM") %>%
  mutate(Station = str_remove(Station, "_TPM$"),
         Station = dplyr::recode(Station, "Z3_2" = "Z4_2", "Z3_V2" = "Z4_V2"),
         TPM = suppressWarnings(as.numeric(TPM)),
         Symbol = str_trim(as.character(Symbol)),
         Pathway = tidy_process(Pathway)) %>%
  filter(Station %in% DESIGN$sqm_id, !is.na(TPM), TPM > 0, !is.na(Symbol)) %>%
  left_join(DESIGN_SQM, by = "Station")
if (!all(DESIGN$sqm_id %in% ncyc$Station))
  warning("NCycDB table: no data for ", paste(setdiff(DESIGN$sqm_id, ncyc$Station), collapse = ", "))

ncyc_sets <- read_excel(CFG$in_ncyc_sets) %>%
  transmute(Pathway = tidy_process(Pathway), Symbol = str_trim(as.character(Symbol))) %>%
  distinct()
unmatched <- setdiff(unique(ncyc$Symbol), ncyc_sets$Symbol)
if (length(unmatched)) warning("Gene families missing from Ncyccompleteness.xlsx: ",
                               paste(unmatched, collapse = ", "))
station_labels <- setNames(as.character(sort(unique(DESIGN$station))), seq_along(unique(DESIGN$station)))

# ---- D2. Gene-family richness (Figure 5A) ------------------------------------
ms("3.5 NCycDB gene families: in the database file / CF / VF (text: 68, 61, 48)",
   c(n_distinct(ncyc_sets$Symbol), n_distinct(ncyc$Symbol[ncyc$Fraction == "Cellular"]),
     n_distinct(ncyc$Symbol[ncyc$Fraction == "Vesicular"])))
n_rich <- ncyc %>% group_by(Station, Label, Fraction, StationNum, x_idx) %>%
  summarise(n_genes = n_distinct(Symbol)) %>% ungroup()
rich_pairs <- PAIRS %>%
  mutate(CF = n_rich$n_genes[match(cf, n_rich$Station)],
         VF = n_rich$n_genes[match(vf, n_rich$Station)],
         VF_pct_of_CF = round(100 * VF / CF, 1)) %>% dplyr::select(station, zone, CF, VF, VF_pct_of_CF)
print(rich_pairs)

p_n_rich <- ggplot(n_rich, aes(x_idx, n_genes, colour = Fraction, group = Fraction)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  scale_x_continuous(breaks = as.integer(names(station_labels)), labels = station_labels) +
  scale_colour_manual(values = PAL$frac_col, guide = "none") +
  labs(x = LAB$stations, y = LAB$n_rich) + theme_tfm()

# ---- D3. Nitrogen compounds (Figure 5B) ---------------------------------------------------
nutrients <- env %>%
  dplyr::select(Station, any_of(c("NO3", "NO2", "NH4", "O.Sat", "Sal"))) %>%
  mutate(across(-Station, ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(Station %in% DESIGN$station) %>% arrange(Station)
print(nutrients)
nut_long <- nutrients %>%
  dplyr::select(Station, any_of(c("NO3", "NO2", "NH4"))) %>%
  pivot_longer(-Station, names_to = "Compound", values_to = "Concentration") %>%
  filter(!is.na(Concentration)) %>%
  group_by(Compound) %>% mutate(Relative = 100 * Concentration / max(Concentration)) %>% ungroup() %>%
  mutate(x_idx = match(Station, sort(unique(DESIGN$station))))
p_nut <- ggplot(nut_long, aes(x_idx, Relative, colour = Compound, group = Compound)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  scale_x_continuous(breaks = as.integer(names(station_labels)), labels = station_labels) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  scale_colour_manual(values = c(NO3 = "#9F00FF", NO2 = "#F8B62D", NH4 = "#5BA664"),
                      labels = c(NO3 = expression(NO[3]^"-"), NO2 = expression(NO[2]^"-"),
                                 NH4 = expression(NH[4]^"+")), name = "Nutrients") +
  labs(x = LAB$stations, y = LAB$n_conc) + theme_tfm()

# ---- D4. Abundance per process (Figure 5C) --------------------------------------------------
# TPM summed per process and scaled to the maximum OF EACH FRACTION
n_abund <- ncyc %>%
  group_by(Fraction, Pathway, x_idx) %>% summarise(TPM_total = sum(TPM)) %>%
  group_by(Fraction, Pathway) %>% mutate(TPM_pct = 100 * TPM_total / max(TPM_total)) %>% ungroup()
p_n_abund <- ggplot(n_abund, aes(x_idx, TPM_pct, colour = Fraction, group = Fraction)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  facet_wrap(~ Pathway, ncol = 2) +
  scale_x_continuous(breaks = as.integer(names(station_labels)), labels = station_labels) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  scale_colour_manual(values = PAL$frac_col, name = "Fraction") +
  labs(x = LAB$stations, y = LAB$n_abund) + theme_tfm()

if (CFG$run_new) {
  n_total <- ncyc %>% group_by(Station) %>% summarise(N_TPM = sum(TPM)) %>%
    left_join(ko_depth, by = "Station") %>% filter(Fraction == "Cellular")
  ct <- suppressWarnings(cor.test(n_total$N_TPM, n_total$Depth, method = "spearman"))
  cat("\n[NEW] Total N-cycle TPM and KO-assigned reads, CF stations:\n")
  print(n_total %>% dplyr::select(Label, N_TPM, Depth) %>% mutate(N_TPM = round(N_TPM, 1)))
  ms("[NEW] Spearman total N-cycle TPM vs depth, CF: rho, P",
     c(round(unname(ct$estimate), 2), signif(ct$p.value, 2)))}

# ---- D5. "Completeness" per process (Figure 5D) -----------------------------------------------
# = percentage of the NCycDB gene families of a process detected in a sample.
n_total_sets <- ncyc_sets %>% dplyr::count(Pathway, name = "n_total")
n_compl <- ncyc %>% distinct(Station, Label, Fraction, Symbol) %>%
  inner_join(ncyc_sets, by = "Symbol", relationship = "many-to-many") %>%
  distinct(Station, Label, Fraction, Pathway, Symbol) %>%
  dplyr::count(Station, Label, Fraction, Pathway, name = "n_detected") %>%
  complete(nesting(Station, Label, Fraction), Pathway = n_total_sets$Pathway,
           fill = list(n_detected = 0)) %>%
  left_join(n_total_sets, by = "Pathway") %>%
  mutate(pct = 100 * n_detected / n_total,
         Pathway = factor(Pathway, levels = rev(sort(unique(Pathway)))),
         Label = factor(str_remove(Label, "^V"), levels = as.character(sort(unique(DESIGN$station)))),
         Fraction_lab = factor(paste(Fraction, "fraction"), levels = paste(FRAC_LEVELS, "fraction")))
save_tab(n_compl, "TableS_nitrogen_gene_families_detected")

print(ncyc %>% inner_join(ncyc_sets %>% filter(str_detect(Pathway, regex("anammox", ignore_case = TRUE))) %>%
                            dplyr::select(Symbol), by = "Symbol") %>%
        group_by(Label) %>% summarise(genes = paste(sort(unique(Symbol)), collapse = ", ")), n = Inf)

p_n_compl <- ggplot(n_compl, aes(Label, Pathway, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  facet_grid(. ~ Fraction_lab, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "#FFF7BC", high = "#B30000", limits = c(0, 100), name = LAB$n_compl) +
  labs(x = LAB$stations, y = NULL) +
  theme_tfm() + theme(axis.ticks = element_blank(), axis.line = element_blank())

# ---- D6. Figure 5  -----------------------------------------------------------
fig5 <- (p_n_rich + p_nut) / p_n_abund / p_n_compl +
  plot_layout(heights = c(1, 4, 1.6)) + plot_annotation(tag_levels = "A") & TAG_THEME
save_fig(fig5, "Figure5", 183, 230)

