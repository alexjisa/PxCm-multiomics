# ============================================================
# Genes detectados con TPM >= 1
# Violin plot: C. melo mock, C. melo inoculated y P. xanthii
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(svglite)
  library(scales)
})

# ------------------------------------------------------------
# Archivos de entrada
# ------------------------------------------------------------

cm_file <- "03_Transcriptomics_Cm_TPM.txt"
px_file <- "03_Transcriptomics_Px_TPM.txt"

# ------------------------------------------------------------
# Leer matrices TPM
# ------------------------------------------------------------

cm <- read_tsv(
  cm_file,
  show_col_types = FALSE
)

px <- read_tsv(
  px_file,
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Columnas de anotación
# ------------------------------------------------------------

meta_cols <- c(
  "Geneid",
  "Chr",
  "Start",
  "End",
  "Strand",
  "Length"
)

# ------------------------------------------------------------
# Identificar columnas de muestras
# ------------------------------------------------------------

cm_sample_cols <- setdiff(
  colnames(cm),
  meta_cols
)

px_sample_cols <- setdiff(
  colnames(px),
  meta_cols
)

# ------------------------------------------------------------
# Convertir columnas TPM a formato numérico
# ------------------------------------------------------------

cm <- cm %>%
  mutate(
    across(
      all_of(cm_sample_cols),
      as.numeric
    )
  )

px <- px %>%
  mutate(
    across(
      all_of(px_sample_cols),
      as.numeric
    )
  )

# ------------------------------------------------------------
# Contar genes con TPM >= 1 en cada muestra de melón
# ------------------------------------------------------------

cm_detected <- cm %>%
  summarise(
    across(
      all_of(cm_sample_cols),
      ~ sum(.x >= 1, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Detected_genes"
  ) %>%
  mutate(
    Group = case_when(
      str_detect(Sample, "_M_") ~ "C. melo Mock",
      str_detect(Sample, "_I_") ~ "C. melo Inoculated",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Group)
  )

# ------------------------------------------------------------
# Contar genes con TPM >= 1 en cada muestra de P. xanthii
# ------------------------------------------------------------

px_detected <- px %>%
  summarise(
    across(
      all_of(px_sample_cols),
      ~ sum(.x >= 1, na.rm = TRUE)
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Sample",
    values_to = "Detected_genes"
  ) %>%
  mutate(
    Group = "P. xanthii"
  )

# ------------------------------------------------------------
# Combinar los tres grupos
# ------------------------------------------------------------

plot_df <- bind_rows(
  cm_detected,
  px_detected
) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "C. melo Mock",
        "C. melo Inoculated",
        "P. xanthii"
      )
    )
  )

# ------------------------------------------------------------
# Resumen estadístico
# ------------------------------------------------------------

summary_df <- plot_df %>%
  group_by(Group) %>%
  summarise(
    Samples = n(),
    Mean = mean(
      Detected_genes,
      na.rm = TRUE
    ),
    Median = median(
      Detected_genes,
      na.rm = TRUE
    ),
    SD = sd(
      Detected_genes,
      na.rm = TRUE
    ),
    Minimum = min(
      Detected_genes,
      na.rm = TRUE
    ),
    Maximum = max(
      Detected_genes,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(summary_df)

write_tsv(
  summary_df,
  "Detected_genes_TPM1_summary.tsv"
)

# ------------------------------------------------------------
# Colores
# ------------------------------------------------------------

group_colors <- c(
  "C. melo Mock" = "#7BC96F",
  "C. melo Inoculated" = "#5DA5DA",
  "P. xanthii" = "#8064A2"
)

# ------------------------------------------------------------
# Figura
# ------------------------------------------------------------

p <- ggplot(
  plot_df,
  aes(
    x = Group,
    y = Detected_genes,
    fill = Group,
    colour = Group
  )
) +
  
  geom_violin(
    trim = FALSE,
    scale = "width",
    width = 0.72,
    linewidth = 0.4,
    alpha = 0.78
  ) +
  
  geom_jitter(
    width = 0.055,
    height = 0,
    size = 1.8,
    alpha = 0.60,
    show.legend = FALSE
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.28,
    linewidth = 0.7,
    colour = "black"
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  scale_colour_manual(
    values = group_colors
  ) +
  
  scale_x_discrete(
    labels = c(
      "C. melo Mock" =
        expression(
          bolditalic("C. melo")~bold("Mock")
        ),
      
      "C. melo Inoculated" =
        expression(
          bolditalic("C. melo")~bold("Inoculated")
        ),
      
      "P. xanthii" =
        expression(
          bolditalic("P. xanthii")
        )
    ),
    expand = expansion(
      mult = c(0.15, 0.15)
    )
  ) +
  
  scale_y_continuous(
    labels = comma,
    expand = expansion(
      mult = c(0.03, 0.06)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Detected genes"
  ) +
  
  guides(
    fill = "none",
    colour = "none"
  ) +
  
  theme_minimal(
    base_size = 10,
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
    
    axis.text.x = element_text(
      face = "bold",
      colour = "black",
      size = 11,
      margin = margin(t = 8)
    ),
    
    axis.text.y = element_text(
      face = "plain",
      colour = "black",
      size = 11
    ),
    
    axis.title.y = element_text(
      face = "bold",
      colour = "black",
      size = 14,
      margin = margin(r = 10)
    ),
    
    axis.ticks = element_blank(),
    
    plot.margin = margin(
      10,
      20,
      10,
      10
    )
  )

# ------------------------------------------------------------
# Mostrar figura
# ------------------------------------------------------------

p

# ------------------------------------------------------------
# Guardar figura
# ------------------------------------------------------------

ggsave(
  filename = "Detected_genes_TPM1_violin.png",
  plot = p,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "Detected_genes_TPM1_violin.pdf",
  plot = p,
  width = 7,
  height = 5,
  units = "in",
  bg = "white"
)

ggsave(
  filename = "Detected_genes_TPM1_violin.svg",
  plot = p,
  width = 7,
  height = 5,
  units = "in",
  bg = "white"
)