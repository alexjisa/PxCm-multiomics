# =============================================================================
# Script:  02_pathway_enrichment.R
# Purpose: Plot KEGG pathway coverage for the top 16 enriched pathways from
#          the MarVis-Pathway output, as a bubble plot (coverage on x, pathway
#          on y, bubble size = putative compound hits, color = enrichment
#          score).
#
# Input:   data/metabolomics/14_Metabolomics_Sets.xlsx
# Output:  results/metabolomics/pathway_enrichment/metabolomics_pathway_coverage.{png,pdf,svg}
# Usage:   Rscript scripts/metabolomics/02_pathway_enrichment.R
#          (run from the project root; see data/README.md to obtain inputs)
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(svglite)
})

# Paths are relative to the project root (see Usage in the header above).
data_dir <- file.path("data", "metabolomics")
results_dir <- file.path("results", "metabolomics", "pathway_enrichment")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(data_dir, "14_Metabolomics_Sets.xlsx")

df <- read_excel(input_file)

exclude_routes <- c(
  "Biosynthesis of secondary metabolites - Cucumis sativus (cucumber)",
  "Metabolic pathways - Cucumis sativus (cucumber)"
)

plot_df <- df %>%
  filter(!`Set name` %in% exclude_routes) %>%
  filter(str_starts(as.character(`Set ID`), "KEGG/")) %>%
  mutate(
    Pathway = `Set name` %>%
      str_remove(" - Cucumis sativus \\(cucumber\\)") %>%
      str_trim(),
    Coverage = 100 * `Entry hits` / Size
  ) %>%
  arrange(desc(Score), desc(`Cpd hits`), desc(Coverage)) %>%
  slice_head(n = 16) %>%
  mutate(
    Pathway = str_wrap(Pathway, width = 42),
    Pathway = factor(Pathway, levels = rev(Pathway))
  )

p <- ggplot(
  plot_df,
  aes(
    x = Coverage,
    y = Pathway,
    size = `Cpd hits`,
    colour = Score
  )
) +
  geom_point(alpha = 0.95) +
  
  scale_colour_gradientn(
    colours = c("blue", "purple", "red"),
    name = "Score",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(0.55, "cm"),
      barheight = unit(3.2, "cm"),
      ticks.colour = "white",
      frame.colour = NA
    )
  ) +
  
  scale_size_continuous(
    name = "Putative metabolites",
    range = c(4, 9),
    breaks = c(5, 10, 15),
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        colour = "black",
        alpha = 0.95
      )
    )
  ) +
  
  scale_y_discrete(
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  
  labs(
    x = "Pathway coverage (%)",
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  
  theme(
    text = element_text(
      family = "Helvetica"
    ),
    
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),
    
    axis.text.y = element_text(
      face = "bold",
      colour = "black",
      size = 8
    ),
    
    axis.text.x = element_text(
      colour = "black",
      size = 8
    ),
    
    axis.title.x = element_text(
      face = "bold",
      colour = "black",
      size = 8
    ),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    legend.box = "vertical",
    
    legend.title = element_text(
      face = "plain",
      size = 8,
      colour = "black"
    ),
    
    legend.text = element_text(
      face = "plain",
      size = 8,
      colour = "black"
    ),
    
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(
  file.path(results_dir, "metabolomics_pathway_coverage.png"),
  p,
  width = 12,
  height = 5,
  dpi = 600
)

ggsave(
  file.path(results_dir, "metabolomics_pathway_coverage.pdf"),
  p,
  width = 12,
  height = 5
)

ggsave(
  file.path(results_dir, "metabolomics_pathway_coverage.svg"),
  p,
  width = 12,
  height = 5
)

p