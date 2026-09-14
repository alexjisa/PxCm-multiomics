suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(svglite)
})

input_file <- "15_Metabolomics_Significant.xlsx"
df <- read_excel(input_file, sheet = "Significant")

sample_cols <- names(df)[str_detect(names(df), "^[0-9]+h_[MI]_R[123]")]

top20 <- df %>%
  mutate(
    DynamicVariance = apply(
      select(., all_of(sample_cols)),
      1,
      var,
      na.rm = TRUE
    )
  ) %>%
  arrange(desc(DynamicVariance)) %>%
  slice_head(n = 20) %>%
  mutate(
    Feature = ifelse(
      is.na(FormerY) | FormerY == "",
      round(`m/z (corrected)`, 3),
      round(as.numeric(FormerY), 3)
    ),
    Feature = make.unique(as.character(Feature))
  )

long_df <- top20 %>%
  select(Feature, all_of(sample_cols)) %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample",
    values_to = "Intensity"
  ) %>%
  mutate(
    Time = as.numeric(str_extract(Sample, "^[0-9]+")),
    Condition = case_when(
      str_detect(Sample, "_M_") ~ "Mock",
      str_detect(Sample, "_I_") ~ "Inoculated"
    )
  )

mean_df <- long_df %>%
  group_by(Feature, Time, Condition) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- mean_df %>%
  group_by(Feature) %>%
  mutate(Zscore = as.numeric(scale(MeanIntensity))) %>%
  ungroup() %>%
  mutate(
    TimeLabel = factor(
      paste0(Time, "h"),
      levels = paste0(sort(unique(Time)), "h")
    )
  )

feature_order <- plot_df %>%
  group_by(Feature) %>%
  summarise(
    MaxAbs = max(abs(Zscore), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(MaxAbs)) %>%
  pull(Feature)

plot_df <- plot_df %>%
  mutate(
    Feature = factor(Feature, levels = rev(feature_order)),
    Condition = factor(Condition, levels = c("Mock", "Inoculated"))
  )

p <- ggplot(plot_df, aes(x = TimeLabel, y = Feature, fill = Zscore)) +
  geom_tile(color = "grey92", linewidth = 0.2) +
  facet_wrap(~ Condition, scales = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Row z-score"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "m/z feature"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "grey40",
      fill = NA,
      linewidth = 0.7
    ),
    panel.spacing = unit(4, "mm"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 9,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12,
      margin = margin(5, 0, 5, 0)
    ),
    legend.title = element_text(
      size = 10,
      face = "plain"
    ),
    legend.text = element_text(size = 9),
    legend.position = "right"
  )

print(p)

ggsave(
  "Top20_dynamic_metabolites_heatmap.png",
  p,
  width = 14,
  height = 5,
  dpi = 600
)

ggsave(
  "Top20_dynamic_metabolites_heatmap.svg",
  p,
  width = 14,
  height = 5
)