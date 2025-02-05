################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="path", help="Input sample file"),
  make_option(c("-g", "--geneFunctions_file"), type="character", default=NULL, metavar="path", help="Gene Functions file."),
  make_option(c("-a", "--annoData_file"), type="character", default=NULL, metavar="path", help="Annotation Data file."),
  make_option(c("-p", "--outprefix"), type="character", default='projectID', metavar="string", help="Output prefix.")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

sampleInput=opt$input_file
geneInput=opt$geneFunctions_file
annoInput=opt$annoData_file
outprefix=opt$outprefix

testing="Y"
if (testing == "Y"){
  sampleInput="sampleData.csv"
  geneInput="geneFunctions.csv"
  annoInput="annoData.csv"
  outprefix="test"
}


if (is.null(sampleInput)){
  print_help(opt_parser)
  stop("Please provide an input file.", call.=FALSE)
}

################################################
################################################
## READ IN FILES##
################################################
################################################
sampleData=read.csv(sampleInput,row.names=1)
annoData=read.csv(annoInput,row.names=1)
geneFunctions=read.csv(geneInput,row.names=1)

################################################
################################################
## Set colors##
################################################
################################################
annoColors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

################################################
################################################
## Create a basic heatmap##

# Normalize the data using z-score normalization
normalizedData <- scale(sampleData)

# Generate basic heatmap
pdf(paste0("basic_heatmap_", outprefix, ".pdf"), width=8, height=6)
pheatmap(
  normalizedData,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D",
  fontsize_row = 6,
  fontsize_col = 8,
  color = colorRampPalette(rev(brewer.pal(n=9, name="RdBu")))(100)
)
dev.off()

################################################
################################################

# Create complex heatmap

# Replace NA, NaN, Inf with zero to prevent errors
normalizedData[is.na(normalizedData)] <- 0
normalizedData[is.nan(normalizedData)] <- 0
normalizedData[is.infinite(normalizedData)] <- 0

low_value <- min(normalizedData, na.rm = TRUE)
high_value <- max(normalizedData, na.rm = TRUE)

# Define the number of colors
num_colors <- 100
color_palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(num_colors)

# Generate a break sequence
# Discrete sequence was not working
breaks_seq <- seq(low_value, high_value, length.out = length(color_palette) + 1)

## MAP GENE FUNCTIONS TO ROW ANNOTATIONS 
gene_function_annotation <- geneFunctions$gene_functions[match(rownames(sampleData), rownames(geneFunctions))]

# Replace NA values with "Unknown"
gene_function_annotation[is.na(gene_function_annotation)] <- "Unknown"

# Convert to factor
row_annotation <- data.frame(gene_functions = factor(gene_function_annotation))
rownames(row_annotation) <- rownames(sampleData)  

# Redefine default colors for gene functions
gene_function_colors <- c(
  "Oxidative_phosphorylation" = "#F46D43",
  "Cell_cycle" = "#708238",
  "Immune_regulation" = "#9E0142",
  "Signal_transduction" = "beige",
  "Transcription" = "violet",
  "Unknown" = "gray"  # Ensure 'Unknown' has a color
)

# Redefine annotation colors
annoColors <- list(
  gene_functions = gene_function_colors,
  Group = c("Disease" = "darkgreen", "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

## GENERATE THE HEATMAP

pdf(paste0("complex_heatmap_", outprefix, ".pdf"), width=10, height=8)

pheatmap(
  normalizedData,  
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D",
  fontsize_row = 6,
  fontsize_col = 8,
  color = color_palette,
  
  # Annotations
  annotation_col = annoData,  
  annotation_row = row_annotation,  
  annotation_colors = annoColors,  
  
  # Hide row names but show column names
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  
  # Set breaks
  breaks = breaks_seq
)

dev.off()

################################################
################################################
