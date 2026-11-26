#!/usr/bin/env Rscript
# Compare original-assembly F1 within median sequence-length groups.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
results <- args[1]
figures <- args[2]
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
lengths <- read.delim(file.path(results, "input_lengths.tsv"))
lengths <- subset(lengths, group == "genome" & scope == "scored" & species != "ALL")
lengths$length_group <- ifelse(lengths$median_bp < 10000, "<10 kb",
                              ifelse(lengths$median_bp <= 100000, "10–100 kb", ">100 kb"))
scores <- read.delim(file.path(results, "metrics.tsv"))
scores <- subset(scores, group == "genome" & stratum == "complete" &
                   tool %in% c("tiberius", "carbon"))
stopifnot(nrow(lengths) == 17, nrow(scores) == 102, all(scores$assessment == "PASS"))
d <- merge(scores, lengths[c("species", "median_bp", "length_group")], by = "species")
paired <- merge(subset(d, tool == "tiberius", select = c(species, level, median_bp, length_group, f1)),
                subset(d, tool == "carbon", select = c(species, level, f1)),
                by = c("species", "level"), suffixes = c("_tiberius", "_carbon"))
paired$delta_pp <- paired$f1_carbon - paired$f1_tiberius
stopifnot(nrow(paired) == 51, !anyDuplicated(paired[c("species", "level")]))
groups <- c("<10 kb", "10–100 kb", ">100 kb")
levels <- c("nucleotide", "exon", "gene")
tests <- list()
for (g in groups) for (lev in levels) {
  x <- subset(paired, length_group == g & level == lev)
  delta <- x$delta_pp
  stopifnot(all(delta != 0), !anyDuplicated(abs(delta)))
  test <- wilcox.test(x$f1_carbon, x$f1_tiberius, paired = TRUE,
                     alternative = "two.sided", exact = TRUE, correct = FALSE)
  # Independently enumerate all signed-rank assignments (at most 2^8).
  ranks <- rank(abs(delta))
  assignments <- expand.grid(rep(list(c(0, 1)), length(delta)))
  possible <- as.vector(as.matrix(assignments) %*% ranks)
  observed <- sum(ranks[delta > 0])
  exact_p <- min(1, 2 * min(mean(possible <= observed), mean(possible >= observed)))
  stopifnot(abs(test$p.value - exact_p) < 1e-12)
  tests[[length(tests) + 1]] <- data.frame(
    length_group = g, level = lev, n_species = nrow(x),
    tiberius_mean_f1 = mean(x$f1_tiberius), carbon_mean_f1 = mean(x$f1_carbon),
    mean_delta_pp = mean(delta), median_delta_pp = median(delta),
    carbon_wins = sum(delta > 0), tiberius_wins = sum(delta < 0),
    V = unname(test$statistic), raw_p = test$p.value)
}
tests <- do.call(rbind, tests)
write.table(paired, file.path(results, "original_assembly_pairs.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(tests, file.path(results, "original_assembly_tests.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
group_n <- tests$n_species[tests$level == "gene"]
group_labels <- paste0(groups, "\nn = ", group_n)
level_labels <- c("Nucleotide", "Exon", "Gene / CDS")
d$length_group <- factor(d$length_group, groups)
d$level <- factor(d$level, levels, level_labels)
d$tool <- factor(d$tool, c("tiberius", "carbon"), c("Tiberius", "Carbon-A"))
d$x <- as.integer(d$length_group) + ifelse(d$tool == "Tiberius", -.12, .12)
tests$x <- match(tests$length_group, groups)
tests$level <- factor(tests$level, levels, level_labels)
tests$label <- sprintf("Mean Δ = %+.1f\nP = %.4f", tests$mean_delta_pp, tests$raw_p)
p <- ggplot(d, aes(x, f1, group = species)) +
  geom_line(colour = "#96999A", linewidth = .3, alpha = .55) +
  geom_point(aes(colour = tool), shape = 16, size = 2, alpha = .88) +
  geom_text(data = tests, aes(x = x, y = 112, label = label), inherit.aes = FALSE,
            size = 2.7, lineheight = 1.05, colour = "#30363B") +
  facet_wrap(~level, nrow = 1) +
  scale_colour_manual(values = c(Tiberius = "#BDA982", `Carbon-A` = "#B76E60")) +
  scale_y_continuous(limits = c(25, 122), breaks = c(25, 50, 75, 100), expand = c(0, 0)) +
  scale_x_continuous(breaks = 1:3, labels = group_labels, limits = c(.52, 3.48), expand = c(0, 0)) +
  labs(title = "Original assemblies", x = NULL, y = "F1 (%)", colour = NULL) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#E8EBED", linewidth = .35),
        panel.spacing = grid::unit(.85, "lines"),
        strip.text = element_text(face = "bold", size = 11, margin = margin(b = 7)),
        legend.position = "top", legend.key.height = grid::unit(.32, "cm"),
        legend.key.width = grid::unit(.55, "cm"),
        legend.margin = margin(t = 0, b = 2),
        plot.title = element_text(face = "bold", size = 15, margin = margin(b = 5)),
        plot.margin = margin(12, 14, 10, 12),
        axis.text = element_text(colour = "#454B50"))
ggsave(file.path(figures, "original_assembly_comparison.png"), p, width = 11.5, height = 4.0, dpi = 220, bg = "white")
ggsave(file.path(figures, "original_assembly_comparison.pdf"), p, width = 11.5, height = 4.0, device = cairo_pdf, bg = "white")
print(tests)
