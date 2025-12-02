# get libraries
library(tidyverse)
library(magrittr)
library(readr)

# VEP output in
vep_data <- read_tsv(
  "scripts/variant_annotation/chr6ch17_mutect2_annotated_with_clinVar.txt",
  comment = "##"
)

# better to set up a theme
theme_course <- function() {
  theme_minimal(base_size = 14) +  # Base size affects all text
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 13, face = "bold"),  # For facets
      legend.position = "bottom"
    )
}

# Split when it contains multiple values
consequences_df <- vep_data %>%
  separate_rows(Consequence, sep="&") %>%
  select(Consequence)

# get chr
vep_data <- vep_data %>%
  separate(Location, into = c("chr", "position"), sep=":") %>% 
  mutate(chr = factor(chr, levels = c("chr6", "chr17")))

# Extract IMPACT from Extra column
vep_data$IMPACT <- gsub(".*IMPACT=([^;]+).*", "\\1", vep_data$Extra)

# get nice data out of the df
extract_field <- function(data, field) {
  gsub(paste0(".*", field, "=([^;]+).*"), "\\1", data$Extra)
}

# combine multiple fields
vep_data <- vep_data %>%
  mutate(
    IMPACT = extract_field(., "IMPACT"),
    SYMBOL = extract_field(., "SYMBOL"),
    VARIANT_CLASS = extract_field(., "VARIANT_CLASS"),
    BIOTYPE = extract_field(., "BIOTYPE")
  )

#  Tip1. Consequence plot, to visualize the distribution of consequences 
# consequence plot
ggplot(vep_data, aes(x=reorder(Consequence, table(Consequence)[Consequence]))) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(title="Distribution of \n Variant Consequences",
       x="Consequence Type",
       y="Count") +
  theme_course()

#  2. Consequence plot, by chromosome 
# consequence plot 
ggplot(vep_data, aes(x=reorder(Consequence, table(Consequence)[Consequence]))) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(title="Distribution of \n Variant Consequences",
       x="Consequence Type",
       y="Count") +
  facet_grid(~chr) +
  theme(
    legend.position = "bottom"
  ) +
  theme_course()

# 3. Can we see the consequence colored by impact 
# Impact + Consequence plot
ggplot(vep_data, aes(x=reorder(Consequence, table(Consequence)[Consequence]), fill=IMPACT)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  labs(title="Variant Consequences \n by Impact",
       x="Consequence Type",
       y="Count") +
  theme(
    legend.position = "bottom"
  )  +
  theme_course()

# 4. Can we see the consequence colored by impact by chromosome 
# Impact + consequence plot x chr
ggplot(vep_data, aes(x=reorder(Consequence, table(Consequence)[Consequence]), fill=IMPACT)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  labs(title="Variant Consequences \n by Impact",
       x="Consequence Type",
       y="Count") +
  facet_grid(~chr) +
  theme(
    legend.position = "bottom"
  )  +
  theme_course()
# 
# Tip5. What are the most common variant classes? 
# variant class plot
vep_data %>%
  count(VARIANT_CLASS) %>%  # Count the classees
  arrange(desc(n)) %>%  # Sort in descending ord
  ggplot(aes(x = reorder(VARIANT_CLASS, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of \n Variant Classes",
       x = "Variant Class",
       y = "Count")  +
  theme_course()
# 6. What are the most common biotypes? 
# biotype distribution
ggplot(vep_data, aes(x=reorder(BIOTYPE, table(BIOTYPE)[BIOTYPE]))) +
  geom_bar(fill="steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title="Distribution of \n Biotypes",
       x="Biotype",
       y="Count") +
  theme_course()

#  Tip7. Bonus: ClinVar, SIFT, ans PolyPhen 
# Extract SIFT
vep_data$sift <- str_extract(vep_data$Extra, "SIFT=.*?;")
vep_data$sift <- str_replace(vep_data$sift, "SIFT=", "") 
vep_data$sift <- str_replace(vep_data$sift, ";", "")
vep_data$sift <- str_replace(vep_data$sift, "\\(.*\\)", "")

# Extract PolyPhen 
vep_data$polyphen <- vep_data$Extra %>%
  str_extract("PolyPhen=.*?;") %>%
  str_replace("PolyPhen=", "") %>%
  str_replace(";", "") %>%
  str_replace("\\(.*\\)", "")

# Extract ClinVar
vep_data$clinvar <- vep_data$Extra %>%
  str_extract("CLIN_SIG=.*?;") %>%
  str_replace("CLIN_SIG=", "") %>%
  str_replace(";", "")

# Distribution of SIFT predictions
ggplot(subset(vep_data, !is.na(sift)), 
       aes(x = reorder(sift, table(sift)[sift]))) +
  geom_bar(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of SIFT Predictions",
       x = "SIFT Prediction",
       y = "Count",
       caption = "Only missense variants are scored by SIFT")  +
  theme_course()

# Distribution of PolyPhen predictions 
ggplot(subset(vep_data, !is.na(polyphen)), 
       aes(x = reorder(polyphen, table(polyphen)[polyphen]))) +
  geom_bar(fill = "darkred") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Distribution of PolyPhen Predictions",
       x = "PolyPhen Prediction",
       y = "Count",
       caption = "Only missense variants are scored by PolyPhen")  +
  theme_course()

# Only for variants with both predictions
p3 <- subset(vep_data, !is.na(sift) & !is.na(polyphen)) %>%
  ggplot(aes(x = sift, fill = IMPACT)) +
  geom_bar() +
  facet_wrap(~polyphen) +
  coord_flip() +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "SIFT and PolyPhen Predictions by Impact",
       x = "SIFT Prediction",
       y = "Count") +
  theme(legend.position = "bottom") +
  theme_course()

# 4. ClinVar significance with impact 
p4 <- subset(vep_data, !is.na(clinvar)) %>%
  ggplot(aes(x = reorder(clinvar, table(clinvar)[clinvar]), fill = IMPACT)) +
  geom_bar() +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "Clinical Significance by Impact",
       x = "ClinVar Clinical Significance",
       y = "Count") +
  theme(legend.position = "bottom") +
  theme_course()

# 5. For missense variants only - comparison of predictions
p5 <- subset(vep_data, Consequence == "missense_variant" & !is.na(sift) & !is.na(polyphen)) %>%
  ggplot(aes(x = sift, fill = polyphen)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "SIFT vs PolyPhen Predictions\nfor Missense Variants",
       x = "SIFT Prediction",
       y = "Count",
       fill = "PolyPhen Prediction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_course()

# Gene-level summary - most frequently mutated genes
p6 <- vep_data %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  count(SYMBOL, chr) %>%
  arrange(desc(n)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(SYMBOL, n), y = n, fill = chr)) +
  geom_col() +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Top 20 Most Frequently Mutated Genes",
       x = "Gene Symbol",
       y = "Number of Variants",
       fill = "Chromosome") +
  theme_course()

# High impact variants by gene
p7 <- vep_data %>%
  filter(IMPACT %in% c("HIGH", "MODERATE") & !is.na(SYMBOL) & SYMBOL != "") %>%
  count(SYMBOL, IMPACT) %>%
  arrange(desc(n)) %>%
  slice_head(n = 30) %>%
  ggplot(aes(x = reorder(SYMBOL, n), y = n, fill = IMPACT)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("HIGH" = "#d62728", "MODERATE" = "#ff7f0e")) +
  labs(title = "Genes with High/Moderate Impact Variants",
       x = "Gene Symbol",
       y = "Number of Variants",
       fill = "Impact Level") +
  theme_course()

print(p3)
print(p4)
print(p5)
print(p6)
print(p7)

# make a summary pf the mutations
summary_table <- vep_data %>%
  summarise(
    total_variants = n(),
    high_impact = sum(IMPACT == "HIGH", na.rm = TRUE),
    moderate_impact = sum(IMPACT == "MODERATE", na.rm = TRUE),
    low_impact = sum(IMPACT == "LOW", na.rm = TRUE),
    with_clinvar = sum(!is.na(clinvar)),
    pathogenic_clinvar = sum(grepl("pathogenic", clinvar, ignore.case = TRUE), na.rm = TRUE),
    with_sift = sum(!is.na(sift)),
    sift_deleterious = sum(grepl("deleterious", sift, ignore.case = TRUE), na.rm = TRUE),
    with_polyphen = sum(!is.na(polyphen)),
    polyphen_damaging = sum(grepl("damaging", polyphen, ignore.case = TRUE), na.rm = TRUE)
  )
