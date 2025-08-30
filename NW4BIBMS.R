#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


# Install if missing:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("pwalign")

library(shiny)
library(Biostrings)
library(pwalign)

# Load the BLOSUM62 matrix from a file
BLOSUM62_df <- read.table("BLOSUM62", header = TRUE, row.names = 1)

# Convert to numeric matrix
# Keep only the 20 standard amino acids
AA <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
BLOSUM62 <- as.matrix(BLOSUM62_df[AA, AA])  # subset rows and columns

# Function for NW alignment and score
nw_align <- function(seq1, seq2) {
  aln <- pairwiseAlignment(seq1, seq2,
                           substitutionMatrix = BLOSUM62,
                           gapOpening = -11, gapExtension = -1,
                           scoreOnly = FALSE)  
  return(aln)
}

# Function to shuffle sequence
shuffle_seq <- function(seq) {
  paste(sample(strsplit(seq, "")[[1]]), collapse = "")
}

# UI
ui <- fluidPage(
  titlePanel("Protein Sequence Alignment with Needleman-Wunsch"),
  # Subtitle below the main title
  tags$p("For BIBMS students – developed by Luis del Peso", 
         style = "font-style: italic; color: #555;"),
  sidebarLayout(
    sidebarPanel(
      textAreaInput("seq1", "Protein sequence 1:", 
                    value = "MTEITAAMVKELRESTGAGMMDCKNALSETQHEWAY"),
      textAreaInput("seq2", "Protein sequence 2:", 
                    value = "MNKMDLVADVAEKTDLSKAKATEVIDAVFA"),
      actionButton("run", "Run Alignment"),
      helpText("Using BLOSUM62 as substitution matrix and -11, -1 as gap open and extension penalties respectively 
(NCBI's defaults for the Needleman-Wunsch Global Align Protein Sequences tool).")
      
    ),
    mainPanel(
      h3("Needleman–Wunsch Alignment Result"),
      verbatimTextOutput("alignment"),
      h4("Alignment Score"),
      verbatimTextOutput("score"),
      h3("Random Shuffle Control Scores (200)"),
      verbatimTextOutput("random_scores"),
      downloadButton("download_scores", "Download Random Scores")
    )
  )
)

# Server
server <- function(input, output) {
  observeEvent(input$run, {
    seq1 <- gsub("[^A-Z]", "", toupper(input$seq1))
    seq2 <- gsub("[^A-Z]", "", toupper(input$seq2))
    
    aln <- nw_align(seq1, seq2)
    
    output$alignment <- renderPrint({
      aln
    })
    
    output$score <- renderPrint({
      score(aln)
    })
    
    random_scores <- replicate(200, {
      s1 <- shuffle_seq(seq1)
      s2 <- shuffle_seq(seq2)
      score(nw_align(s1, s2))
    })
    
    output$random_scores <- renderPrint({
      random_scores
    })
    
    output$download_scores <- downloadHandler(
      filename = function() {
        paste0("random_scores_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Convert the vector to a data frame with one column
        write.csv(data.frame(Score = random_scores), file, row.names = FALSE)
      })
      
  })
}

# Run app
shinyApp(ui = ui, server = server)
