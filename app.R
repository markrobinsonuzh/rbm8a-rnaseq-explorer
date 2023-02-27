suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(ggrepel)
  library(GGally)
  library(dplyr)
  library(tidyr)
  library(Gviz)
  library(ComplexHeatmap)
  library(reshape2)
  library(circlize)
  library(gridExtra)
})

topdir <- "./"

## Load help functions and summarized data
source(file.path(topdir, "plot_tracks_simple.R"))
res <- readRDS(file.path(topdir, "shiny_results.rds"))
options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)
options(ucscChromosomeNames=FALSE)


col_fun_dex <- colorRamp2(c(0, .5, 1), 
                          c("blue", "white", "red"))


# hack to remove the d3d5 sample (lost line)
# ----------------------------------------------
res$metadata <- res$metadata %>% dplyr::filter(group != "24h_rbm8a_d3d5")

res$pca <- res$pca %>% dplyr::filter(group != "24h_rbm8a_d3d5")

res$cpms <- res$cpms[,-grep("20170530.A-8", colnames(res$cpms))]

del <- grep("20170530.A-8_trimmed_Aligned.sortedByCoord.out", names(res$bw_files))
res$bw_files <- res$bw_files[-del]

res$res_long <- res$res_long %>% dplyr::filter(!grepl("d3d5", contrast))

del <- grep("d3d5", colnames(res$res_wide))
res$res_wide <- res$res_wide[,-del]

del <- grep("d3d5", colnames(res$camerares_wide))
res$camerares_wide <- res$camerares_wide[,-del]

del <- grep("d3d5", colnames(res$dexseqres_wide))
res$dexseqres_wide <- res$dexseqres_wide[,-del]

del <- grep("d3d5", colnames(res$suppares_wide))
res$suppares_wide <- res$suppares_wide[,-del]

del <- grep("d3d5", colnames(res$dexseqdtures_wide))
res$dexseqdtures_wide <- res$dexseqdtures_wide[,-del]
# ----------------------------------------------
# ----------------------------------------------

res$dex_w <- res$dexseqres_wide %>%
  group_by(gene_id) %>%
  dplyr::summarize(symbol = symbol[1],
                   rbm8ad5mm.vs.RGB.24h = as.integer(any(padj_gene.rbm8ad5mm.vs.RGB.24h < .05)),
                   #rbm8ad3d5.vs.RGB.24h = as.integer(any(padj_gene.rbm8ad3d5.vs.RGB.24h < .05)),
                   #rbm8ad5mm.vs.rbm8ad3d5.24h = as.integer(any(padj_gene.rbm8ad5mm.vs.rbm8ad3d5.24h < .05)),
                   rbm8ad5mm.vs.RGB.budstage = as.integer(any(padj_gene.rbm8ad5mm.vs.RGB.budstage < .05))) %>%
  dplyr::select(gene_id, symbol, rbm8ad5mm.vs.RGB.24h,
                #rbm8ad3d5.vs.RGB.24h, rbm8ad5mm.vs.rbm8ad3d5.24h,
                rbm8ad5mm.vs.RGB.budstage)

res$suppa_w <- res$suppares_wide %>%
  group_by(gene, etype) %>%
  dplyr::summarize(symbol = symbol[1],
                   rbm8ad5mm.vs.RGB.24h = as.integer(any(padj.rbm8ad5mm.vs.RGB.24h < .05)),
                   #rbm8ad3d5.vs.RGB.24h = as.integer(any(padj.rbm8ad3d5.vs.RGB.24h < .05)),
                   #rbm8ad5mm.vs.rbm8ad3d5.24h = as.integer(any(padj.rbm8ad5mm.vs.rbm8ad3d5.24h < .05)),
                   rbm8ad5mm.vs.RGB.budstage = as.integer(any(padj.rbm8ad5mm.vs.RGB.budstage < .05))) %>%
  dplyr::select(gene, symbol, etype, rbm8ad5mm.vs.RGB.24h,
                rbm8ad5mm.vs.RGB.budstage)
  #dplyr::select(gene, symbol, etype, rbm8ad5mm.vs.RGB.24h,
  #              rbm8ad3d5.vs.RGB.24h, rbm8ad5mm.vs.rbm8ad3d5.24h,
  #              rbm8ad5mm.vs.RGB.budstage)


res$dex_labels <- dex_w <- res$dexseqres_wide %>%
  group_by(gene_id) %>%
  summarize(gene = paste0("gene: (",
                          as.integer(any(padj_gene.rbm8ad5mm.vs.RGB.24h < .05)), ",",
                          #as.integer(any(padj_gene.rbm8ad3d5.vs.RGB.24h < .05)), ",",
                          #as.integer(any(padj_gene.rbm8ad5mm.vs.rbm8ad3d5.24h < .05)), ",",
                          as.integer(any(padj_gene.rbm8ad5mm.vs.RGB.budstage < .05)), ")"),
            symbol = symbol[1],
            intron = paste0("intron: (",
                            sum(padj.rbm8ad5mm.vs.RGB.24h[ftype=="intron"] < .05), ",",
                            #sum(padj.rbm8ad3d5.vs.RGB.24h < .05), ",",
                            #sum(padj.rbm8ad5mm.vs.rbm8ad3d5.24h < .05), ",",
                            sum(padj.rbm8ad5mm.vs.RGB.budstage < .05), ")"),
            exon = paste0("exon: (",
                          sum(padj.rbm8ad5mm.vs.RGB.24h < .05), ",",
                          #sum(padj.rbm8ad3d5.vs.RGB.24h < .05), ",",
                          #sum(padj.rbm8ad5mm.vs.rbm8ad3d5.24h < .05), ",",
                          sum(padj.rbm8ad5mm.vs.RGB.budstage < .05), ")"), 
            min_p_gene = min(padj_gene.rbm8ad5mm.vs.RGB.24h,
                             #padj_gene.rbm8ad3d5.vs.RGB.24h,
                             #padj_gene.rbm8ad5mm.vs.rbm8ad3d5.24h,
                             padj_gene.rbm8ad5mm.vs.RGB.budstage, na.rm = TRUE),
            min_p_feature = min(padj.rbm8ad5mm.vs.RGB.24h,
                             #   padj.rbm8ad3d5.vs.RGB.24h,
                             #   padj.rbm8ad5mm.vs.rbm8ad3d5.24h,
                                padj.rbm8ad5mm.vs.RGB.budstage, na.rm = TRUE)) %>%
  mutate(label=paste0(symbol, " (", gene_id, ")\n",
         "(d5.vs.RGB.24h, d5mm.vs.RGB.budstage)\n",
         #"(d5.vs.RGB.24h, d3d5.vs.RGB.24h, d5mm.vs.d3d5.24h, d5mm.vs.RGB.budstage)\n",
         gene, " / ", intron, " / ", exon,"\n",
         "min_P (gene: ", min_p_gene, "; feature: ", min_p_feature, ")\n")) %>%
  select(gene_id, symbol, label)

comps <- setdiff(colnames(res$dex_w), c("gene_id", "symbol"))
stopifnot( all(comps %in% colnames(res$suppa_w)) )


## Define app
summary_app <- function(res) {
  options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)
  
  p_layout <- function(request) {
    shinydashboard::dashboardPage (
      skin = "blue", 
      
      shinydashboard::dashboardHeader(title = "zebrafish RNA-seq, p2452", titleWidth = 300),
      
      shinydashboard::dashboardSidebar(
        textInput(inputId = "sel.gene", label = "Selected gene"),
        bookmarkButton()
      ),
      
      shinydashboard::dashboardBody(fluidRow(
        shinydashboard::tabBox(
          width = 12, 
          
          tabPanel("About",
                   includeMarkdown(paste0(topdir, "shiny/about_app.md"))),
          
          tabPanel("PCA",
                   uiOutput("pca.plot.ui"),
                   selectInput("pcx", "x-axis", choices = paste0("PC", 1:7),
                               selected = "PC1"),
                   selectInput("pcy", "y-axis", choices = paste0("PC", 1:7),
                               selected = "PC2")),
          
          tabPanel("Gene CPM plot",
                   div(style = "position:relative",
                       uiOutput("cpm.plot.ui"),
                       uiOutput("cpm.hover.info"))),
          
          tabPanel("Coverage plot",
                   uiOutput("coverage.plot.ui")),
          
          tabPanel("Genome-wide DGE plot",
                   fluidRow(
                     column(5, textInput(inputId = "sel.chrom", label = "Chromosome")),
                     column(5, div(style = "margin-top: 25px;", textOutput("chr.length.text")))
                     #column(5, textOutput("chr.length.text"))
                   ),
                   #tags$style(type = 'text/css', "#chr.length.text { width:100%; margin-top: 25px;}"),
                   fluidRow(
                     column(4, textInput(inputId = "sel.start", label = "Start position")),
                     column(4, textInput(inputId = "sel.end", label = "End position")),
                     column(4, numericInput(inputId = "fdrbaseline", 
                                            label = "FDR threshold", value = 0.05, 
                                            min = 0, max = 1, step = 0.001))
                   ),
                   uiOutput("genomedge.plot.ui")),
          
          tabPanel("DGE volcano plots",
                   div(style = "position:relative",
                       uiOutput("volcano.plot.ui"),
                       uiOutput("volcano.hover.info"))),
          
          tabPanel("DGE result table",
                   DT::dataTableOutput("result.table.edger"),
                   p(class = "text-center", 
                     downloadButton("result.table.edger.download", "Download"))),
          
          tabPanel("DGE GO camera result table",
                   DT::dataTableOutput("result.table.camera"),
                   p(class = "text-center",
                     downloadButton("result.table.camera.download", "Download"))),
          
          tabPanel("Gene set heatmap",
                   selectInput(inputId = "sel.geneset.samples",
                               label = "Included samples", 
                               choices = colnames(res$cpms), multiple = TRUE,
                               width = "100%", selected = colnames(res$cpms)),
                   selectInput(inputId = "sel.geneset", label = "Selected gene set",
                               choices = res$camerares_wide$geneset, width = "100%"),
                   sliderInput("height", "plot height", min = 500, max = 3000, value = 700),
                   sliderInput("width", "plot width", min = 700, max = 1000, value = 750),
                   uiOutput("geneset.heatmap.ui")),
          
          tabPanel("DEXSeq summary heatmap",
                   selectInput(inputId = "sel.comparisons.dexseq",
                               label = "Included comparisons", 
                               choices = comps, multiple = TRUE,
                               width = "100%", selected = comps),
                   selectInput(inputId = "sel.geneset.dexseq", label = "Selected gene set",
                               choices = res$camerares_wide$geneset, width = "100%"),
                   sliderInput("height.dexseq", "plot height", min = 500, max = 3000, value = 700),
                   sliderInput("width.dexseq", "plot width", min = 700, max = 1000, value = 750),
                   uiOutput("dexseq.heatmap.ui")),

          tabPanel("DEXSeq result table",
                   DT::dataTableOutput("result.table.dexseq"),
                   p(class = "text-center", 
                     downloadButton("result.table.dexseq.download", "Download"))),

          tabPanel("DEXSeq bin-count plot",
                   selectInput(inputId = "sel.dexgroups",
                               label = "Included groups", 
                               choices = unique(res$metadata$group), multiple = TRUE,
                               width = "100%", selected = unique(res$metadata$group)),
                   uiOutput("dexseq.binplot.ui")),

          tabPanel("DEXSeq DTU bin-count plot",
                   selectInput(inputId = "sel.dexdtugroups",
                               label = "Included groups", 
                               choices = unique(res$metadata$group), multiple = TRUE,
                               width = "100%", selected = unique(res$metadata$group)),
                   uiOutput("dexseqdtu.binplot.ui")),
          
          tabPanel("SUPPA summary heatmap",
                   selectInput(inputId = "sel.comparisons.suppa",
                               label = "Included comparisons", 
                               choices = comps, multiple = TRUE,
                               width = "100%", selected = comps),
                   selectInput(inputId = "sel.geneset.suppa", label = "Selected gene set",
                               choices = res$camerares_wide$geneset, width = "100%"),
                   sliderInput("height.suppa", "plot height", min = 500, max = 3000, value = 700),
                   sliderInput("width.suppa", "plot width", min = 700, max = 1000, value = 750),
                   uiOutput("suppa.heatmap.ui")),
          
          tabPanel("SUPPA result table",
                   DT::dataTableOutput("result.table.suppa"),
                   p(class = "text-center", 
                     downloadButton("result.table.suppa.download", "Download"))),
          
          tabPanel("DEXSeq DTU result table",
                   DT::dataTableOutput("result.table.dexseq.dtu"),
                   p(class = "text-center", 
                     downloadButton("result.table.dexseq.dtu.download", "Download")))
        )
      ))
    )
  }

  server_function <- function(input, output, session) {
    options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)
    
    output$chr.length.text <- 
      renderText(paste0("The length of this chromosome is ", 
                        res$chromlengths$length[which(res$chromlengths$chromosome == 
                                                        input$sel.chrom)], " bp"))
    
    # ================== Result table camera ======================
    output$result.table.camera <- DT::renderDataTable({
      res$camerares_wide
    },
    filter = "top",
    rownames = FALSE, 
    options = list(scrollX = TRUE))
    
    output$result.table.camera.download <- 
      downloadHandler("camera_results_filtered.csv", content = function(file) {
        s.camera = input$result.table.camera_rows_all
        write.csv(res$camerares_wide[s.camera, , drop = FALSE], file)
      })
    
    # ================== Result table edgeR =======================
    output$result.table.edger <- DT::renderDataTable({
      res$res_wide %>%
        dplyr::mutate(gene_biotype = factor(gene_biotype))
    },
    filter = "top",
    rownames = FALSE, 
    options = list(scrollX = TRUE))
    
    output$result.table.edger.download <- 
      downloadHandler("DGE_results_filtered.csv", content = function(file) {
      s.edger = input$result.table.edger_rows_all
      write.csv(res$res_wide[s.edger, , drop = FALSE], file)
    })
    
    # ================== Result table DEXSeq ======================
    output$result.table.dexseq <- DT::renderDataTable({
      res$dexseqres_wide %>%
        dplyr::mutate(gene_biotype = factor(gene_biotype),
                      ftype = factor(ftype),
                      seqnames = factor(seqnames),
                      strand = factor(strand),
                      U12spliceintron = factor(U12spliceintron),
                      trueintron = factor(trueintron),
                      trueexon = factor(trueexon))
    },
    filter = "top",
    rownames = FALSE, 
    options = list(scrollX = TRUE))
 
    output$result.table.dexseq.download <- 
      downloadHandler("DEXSeq_results_filtered.csv", content = function(file) {
        s.dexseq = input$result.table.dexseq_rows_all
        write.csv(res$dexseqres_wide[s.dexseq, , drop = FALSE], file)
      })
    
    # ================= Result table DEXSeq DTU ====================
    output$result.table.dexseq.dtu <- DT::renderDataTable({
      res$dexseqdtures_wide %>%
        dplyr::mutate(gene_biotype = factor(gene_biotype))
    },
    filter = "top",
    rownames = FALSE, 
    options = list(scrollX = TRUE))
    
    output$result.table.dexseq.dtu.download <- 
      downloadHandler("DEXSeq_DTU_results_filtered.csv", content = function(file) {
        s.dexseq.dtu = input$result.table.dexseq.dtu_rows_all
        write.csv(res$dexseqdtures_wide[s.dexseq.dtu, , drop = FALSE], file)
      })
    
    # ================== Result table SUPPA =======================
    output$result.table.suppa <- DT::renderDataTable({
      res$suppares_wide %>%
        dplyr::mutate(gene_biotype = factor(gene_biotype),
                      etype = factor(etype))
    },
    filter = "top",
    rownames = FALSE, 
    options = list(scrollX = TRUE))
    
    output$result.table.suppa.download <- 
      downloadHandler("SUPPA_results_filtered.csv", content = function(file) {
        s.suppa = input$result.table.suppa_rows_all
        write.csv(res$suppares_wide[s.suppa, , drop = FALSE], file)
      })
    
    # ===================== Coverage plot =========================
    output$coverage.plot <- renderPlot({
      if (input$sel.gene == "")
        return(NULL)
      else {
        idshow <- res$dex_labels %>% 
          dplyr::filter(gene_id == input$sel.gene | symbol == input$sel.gene) %>%
          pull(label)
        if(length(idshow)==0)
          idshow <- NULL
        plot_tracks(mygene = input$sel.gene, genemodels = res$gene_models, 
                    genemodels2 = res$genemodelsexcl, 
                    gtf_file = NULL, rnaseq_datafiles = res$bw_files, 
                    rnaseq_condition = res$condition, show_chr = NULL, 
                    min_coord = NULL, max_coord = NULL,
                    cex.main = 1, idshow = idshow,
                    pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
        
      }
    })
    
    output$coverage.plot.ui <- renderUI({
      plotOutput("coverage.plot", height = "800px")
    })
    
    # =================== Genome-wide DGE plot ===============
    output$genomedge.plot <- renderPlot({
      if (input$sel.chrom == "")
        return(NULL)
      else {
        if (input$sel.start == "") startpos <- 1
        else startpos <- as.numeric(as.character(input$sel.start))
        if (input$sel.end == "") 
          endpos <- as.numeric(as.character(res$chromlengths$length[
            which(res$chromlengths$chromosome == input$sel.chrom)]))
        else endpos <- as.numeric(as.character(input$sel.end))
        if (endpos - startpos <= 500000) doshow <- TRUE
        else doshow <- FALSE
        
        nms <- grep("mlog10FDR", colnames(mcols(res$gtfgene)), value = TRUE)
        fdrTracks <- lapply(1:length(nms), function(i) {
          ir1 <- res$gtfgene
          strand(ir1) <- "*"
          mcols(ir1) <- mcols(ir1)[, nms[i], drop = FALSE]
          grps <- gsub("mlog10FDR\\.", "", colnames(mcols(ir1)))
          
          assign(paste0("fdrtr", i), 
                 DataTrack(range = ir1,
                           type = "histogram",
                           name = gsub("mlog10FDR\\.", "", nms[i]),
                           groups = grps, 
                           alpha = 1, fill = i,
                           ylim = c(-4.5, 4.5),
                           baseline = c(-log10(input$fdrbaseline), 0, 
                                        log10(input$fdrbaseline))))
        })
        fdrplotTracks <- c(fdrTracks, list(GenomeAxisTrack(), res$genetrack))
        
        plotTracks(fdrplotTracks, chromosome = input$sel.chrom,
                   showOverplotting = TRUE, 
                   from = startpos,
                   to = endpos, showId = doshow,
                   col.title = "black", 
                   main = "-log10(FDR) * sign(logFC)")
      }
    })
    
    output$genomedge.plot.ui <- renderUI({
      plotOutput("genomedge.plot", height = "800px")
    })
    
    # ===================== PCA plot =========================
    output$pca.plot <- renderPlot({
      if (!is.null(input$pcx) & !is.null(input$pcy)) {
        ggplot(res$pca, aes_string(x = input$pcx, y = input$pcy, color = "trttime", 
                                   label = "sample_id")) +
          geom_point(size = 7) + geom_label_repel() + theme_bw() + 
          xlab(paste0(input$pcx, res$pcavar[input$pcx])) + 
          ylab(paste0(input$pcy, res$pcavar[input$pcy])) + 
          coord_fixed() + theme(axis.text = element_text(size = 14),
                                axis.title = element_text(size = 16)) + 
          scale_colour_manual(values = c(`rbm8a_d3d5.24h` = "#DC050C", 
                                         `rbm8a_d5--.24h` = "#F1932D", 
                                         `rbm8a_d5--.budstage` = "#B17BA6", 
                                         `RGB.24h` = "#7BAFDE", `RGB.budstage` = "black"))
          #scale_colour_manual(values = c(`rbm8a_d3d5.24h` = "#DC050C",
          #                               `rbm8a_d5--.24h` = "#7BAFDE",
          #                               `rbm8a_d5--.budstage` = "#B17BA6",
          #                               `RGB.24h` = "#F1932D", `RGB.budstage` = "black"))

      } else {
        NULL
      }
    })

    output$pca.plot.ui <- renderUI({
      plotOutput("pca.plot", height = "600px")
    })
    
    # =================== Gene set heatmap =======================
    output$geneset.heatmap <- renderPlot(
      width = function() input$width,
      height = function() input$height,
      {
      if (!is.null(input$sel.geneset) && input$sel.geneset %in% 
          res$camerares_wide$geneset) {
        genes_in_set <- res$gene_info$gene[res$gene_info$symbol %in% 
                                             strsplit(res$camerares_wide$genes[
                                               res$camerares_wide$geneset == 
                                                 input$sel.geneset
                                             ], " ")[[1]]]
        tmp <- res$cpms[as.character(genes_in_set), 
                        as.character(input$sel.geneset.samples), 
                        drop = FALSE]
        rownames(tmp) <- res$gene_info$symbol[match(rownames(tmp), 
                                                    res$gene_info$gene)]
        colAnnot <- ComplexHeatmap::HeatmapAnnotation(
          group = res$metadata$trttime[match(colnames(tmp), res$metadata$sample_id)], 
          col = list(group = c(`rbm8a_d3d5.24h` = "#DC050C", 
                               `rbm8a_d5--.24h` = "#F1932D", 
                               `rbm8a_d5--.budstage` = "#B17BA6", 
                               `RGB.24h` = "#7BAFDE", `RGB.budstage` = "black"))
          #col = list(group = c(`rbm8a_d3d5.24h` = "#DC050C", 
          #                     `rbm8a_d5--.24h` = "#7BAFDE", 
          #                     `rbm8a_d5--.budstage` = "#B17BA6", 
          #                     `RGB.24h` = "#F1932D", `RGB.budstage` = "black"))
        )
        
        ComplexHeatmap::Heatmap(t(scale(t(log2(tmp + 1)), center = TRUE, 
                                          scale = FALSE)),
                                name = "log2(CPM + 1)\ncentered",
                                bottom_annotation = colAnnot)
      } else {
        NULL
      }
    })
    
    output$geneset.heatmap.ui <- renderUI({
      plotOutput("geneset.heatmap")
    })
    
    # =================== DEXSeq heatmap =======================
    output$dexseq.heatmap <- renderPlot(
      width = function() input$width.dexseq,
      height = function() input$height.dexseq,
      {
        if (!is.null(input$sel.geneset.dexseq) && input$sel.geneset.dexseq %in% 
            res$camerares_wide$geneset) {
          
          gs <- res$camerares_wide$geneset == input$sel.geneset.dexseq
          genes_in_set <- strsplit(res$camerares_wide$genes[gs], " ")[[1]]
          
          tmp <- res$dex_w %>% 
            dplyr::filter(symbol %in% genes_in_set)
          rn <- tmp$symbol
          if(is.null(input$sel.comparisons.dexseq)) {
            return(NULL)
          }
          tmp <- tmp %>% 
            dplyr::select(-gene_id,-symbol) %>% 
            as.matrix()
          tmp <- tmp[,input$sel.comparisons.dexseq,drop=FALSE]
          rownames(tmp) <- rn
          # tmp <- tmp + rnorm(prod(dim(tmp)), sd = .001)
          o <- order(rowSums(tmp, na.rm=TRUE))
          ComplexHeatmap::Heatmap(tmp[o,], name = "DEXSeq change",
                                  col = col_fun_dex,
                                  cluster_columns = FALSE, 
                                  cluster_rows = FALSE, #column_names_rot = 60,
                                  show_heatmap_legend = FALSE,
                                  column_names_max_height = unit(90,"mm"),
                                  column_title = paste0(input$sel.geneset.dexseq, "\n",
                                                     "(red=change, blue=no-change, grey=NA)"))
        } else {
          NULL
        }
      })
    
    output$dexseq.heatmap.ui <- renderUI({
      plotOutput("dexseq.heatmap")
    })   
    

    # =================== SUPPA heatmap =======================
    output$suppa.heatmap <- renderPlot(
      width = function() input$width.suppa,
      height = function() input$height.suppa,
      {
        if (!is.null(input$sel.geneset.suppa) && input$sel.geneset.suppa %in% 
            res$camerares_wide$geneset) {
          
          gs <- res$camerares_wide$geneset == input$sel.geneset.suppa
          genes_in_set <- strsplit(res$camerares_wide$genes[gs], " ")[[1]]
          
          tmp <- res$suppa_w %>% 
            dplyr::filter(symbol %in% genes_in_set)
          etype <- tmp$etype
          symbol <- tmp$symbol
          rn <- paste0(etype,": ",symbol)

          if(is.null(input$sel.comparisons.suppa)) {
            return(NULL)
          }
          tmp <- tmp %>% ungroup %>%
            dplyr::select(-gene,-symbol, -etype) %>%
            as.matrix()
          tmp <- tmp[,input$sel.comparisons.suppa,drop=FALSE]
          rownames(tmp) <- rn
          # tmp <- tmp + rnorm(prod(dim(tmp)), sd = .001)
          o <- order(rowSums(tmp, na.rm=TRUE), etype, symbol)
          ComplexHeatmap::Heatmap(tmp[o,], name = "SUPPA change",
                                  col = col_fun_dex,
                                  cluster_columns = FALSE, 
                                  cluster_rows = FALSE, #column_names_rot = 60,
                                  show_heatmap_legend = FALSE,
                                  column_names_max_height = unit(90,"mm"),
                                  column_title = paste0(input$sel.geneset.suppa, "\n",
                                                        "(red=change, blue=no-change, grey=NA)"))
        } else {
          NULL
        }
      })
    
    output$suppa.heatmap.ui <- renderUI({
      plotOutput("suppa.heatmap")
    })   
    
    # ===================== DEXseq bin-count plot =========================
    output$dexseq.binplot <- renderPlot({
      if (input$sel.gene == "")
        return(NULL)
      else {
        id <- res$dex_labels %>% 
          dplyr::filter(gene_id == input$sel.gene | symbol == input$sel.gene) %>%
          pull(gene_id)
        make_simple_plot(id, subset_cols = input$sel.dexgroups)
      }
    })
    output$dexseq.binplot.ui <- renderUI({
      plotOutput("dexseq.binplot", height = "800px")
    })   
    

    # ===================== DEXseq DTU bin-count plot =========================
    output$dexseqdtu.binplot <- renderPlot({
      if (input$sel.gene == "")
        return(NULL)
      else {
        id <- res$dex_labels %>% 
          dplyr::filter(gene_id == input$sel.gene | symbol == input$sel.gene) %>%
          pull(gene_id)
        make_proportions_plot(id, subset_cols = input$sel.dexdtugroups)
      }
    })
    output$dexseqdtu.binplot.ui <- renderUI({
      plotOutput("dexseqdtu.binplot", height = "800px")
    })   

        
    # ===================== Volcano plot =========================
    output$volcano.plot.ui <- renderUI({
      plotOutput("volcano.plot", hover = "volcano_hover",
                 click = "volcano_click", 
                 width = "90%", height = "500px")
    })

    output$volcano.plot <- renderPlot({
      p <- ggplot(res$res_long, aes(x = logFC, y = mlog10PValue, 
                                    col = as.character(FDR <= 0.05))) +
        theme_bw() + ylab("-log10(PValue)") + facet_grid(~contrast) + 
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                           name = "FDR <= 0.05")
      if (input$sel.gene == "") return(p + geom_point(size = 1))
      else {
        if (!(tolower(input$sel.gene) %in% c(gsub("\\.[0-9]{1,2}$", "", tolower(res$res_long$gene)),
                                             tolower(res$res_long$symbol))))
          return(p + geom_point(size = 1))
        id <- (res$res_long %>%
                 dplyr::filter(gsub("\\.[0-9]{1,2}$", "", tolower(gene)) == tolower(input$sel.gene) |
                                 tolower(symbol) == tolower(input$sel.gene)))$gene
        if (is.null(id)) return(p + geom_point(size = 1))
        p + geom_point(size = 1, alpha = 0.1) +
          geom_point(data = res$res_long %>% dplyr::filter(gene %in% id), size = 4, pch = 21,
                     col = "yellow", aes(fill = as.character(FDR <= 0.05))) +
          guides(fill = FALSE) + 
          scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black"), 
                            name = "FDR <= 0.05")
      }
    })
    
    ## Clicking in the volcano plot modifies the selected gene
    pp <- reactive(nearPoints(res$res_long, input$volcano_click, threshold = 5, 
                              maxpoints = 1, panelvar1 = "contrast"))
    observeEvent(input$volcano_click, {
      updateTextInput(session, "sel.gene", value = gsub("\\.[0-9]{1,2}$", "", pp()$gene))
    })
    
    output$volcano.hover.info <- renderUI({
      volcano_hover <- input$volcano_hover
      res1 <- nearPoints(res$res_long, volcano_hover, threshold = 5,
                         maxpoints = 1, panelvar1 = "contrast")
      left_pct <- (volcano_hover$x - volcano_hover$domain$left)/
        (volcano_hover$domain$right - volcano_hover$domain$left)
      top_pct <- (volcano_hover$domain$top - volcano_hover$y)/
        (volcano_hover$domain$top - volcano_hover$domain$bottom)
      left_px <- volcano_hover$range$left +
        left_pct * (volcano_hover$range$right - volcano_hover$range$left)
      top_px <- volcano_hover$range$top +
        top_pct * (volcano_hover$range$bottom - volcano_hover$range$top)
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      wellPanel(
        style = style,
        p(HTML(paste0("<b> ID: </b>", res1$gene, "<br/>",
                      "<b> Gene: </b>", res1$symbol, "<br/>",
                      "<b> logFC: </b>", round(res1$logFC, 3), "<br/>",
                      "<b> PValue: </b>", signif(res1$PValue, 3), "<br/>",
                      "<b> FDR: </b>", signif(res1$FDR, 3), "<br/>")))
      )
    })

    # ===================== CPM plot =========================
    output$cpm.plot.ui <- renderUI({
      plotOutput("cpm.plot", hover = "cpm_hover", width = "90%", height = "500px")
    })

    output$cpm.plot <- renderPlot({
      if (input$sel.gene == "") return(NULL)
      else {
        if (!(tolower(input$sel.gene) %in% c(gsub("\\.[0-9]{1,2}$", "", tolower(res$gene_info$gene)),
                                             tolower(res$gene_info$symbol)))) return(NULL)
        id <- (res$gene_info %>% dplyr::filter(gsub("\\.[0-9]{1,2}$", "", tolower(gene)) == tolower(input$sel.gene) |
                                          tolower(symbol) == tolower(input$sel.gene)))$gene
        if (is.null(id)) return(NULL)
        else {
          df <- reshape2::melt(data.frame(t(res$cpms[match(id, rownames(res$cpms)), , drop = FALSE]),
                                          sample = colnames(res$cpms),
                                          group = res$metadata[match(colnames(res$cpms),
                                                                     res$metadata$sample_id), "trttime"]))
          ggplot(df, aes(x = sample, y = value, group = variable, col = group)) +
            geom_line(col = "black") + geom_point(size = 3) +
            theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
                               legend.position = "bottom",
                               axis.text.y = element_text(size = 14),
                               axis.title.y = element_text(size = 16)) +
            guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + xlab("") + 
            ylab("CPM") + 
            scale_colour_manual(values = c(`rbm8a_d3d5.24h` = "#DC050C", 
                                         `rbm8a_d5--.24h` = "#F1932D",       
                                         `rbm8a_d5--.budstage` = "#B17BA6",
                                         `RGB.24h` = "#7BAFDE", `RGB.budstage` = "black"))
          #scale_colour_manual(values = c(`rbm8a_d3d5.24h` = "#DC050C",
          #                               `rbm8a_d5--.24h` = "#7BAFDE",
          #                               `rbm8a_d5--.budstage` = "#B17BA6",
          #                               `RGB.24h` = "#F1932D", `RGB.budstage` = "black"))
        }
      }
    })
    
    output$cpm.hover.info <- renderUI({
      if (input$sel.gene == "") return(NULL)
      else {
        if (!(tolower(input$sel.gene) %in% c(gsub("\\.[0-9]{1,2}$", "", tolower(res$gene_info$gene)),
                                             tolower(res$gene_info$symbol)))) return(NULL)
        id <- (res$gene_info %>% dplyr::filter(gsub("\\.[0-9]{1,2}$", "", tolower(gene)) == tolower(input$sel.gene) |
                                      tolower(symbol) == tolower(input$sel.gene)))$gene
        if (is.null(id)) return(NULL)
        else {
          cpm_hover <- input$cpm_hover
          df <- reshape2::melt(data.frame(t(res$cpms[match(id, rownames(res$cpms)), , drop = FALSE]),
                                          sample = colnames(res$cpms),
                                          group = res$metadata[match(colnames(res$cpms), 
                                                                     res$metadata$sample_id), "trttime"]))
          res2 <- nearPoints(df, cpm_hover, threshold = 5, maxpoints = 1)
          res2$genename <- res$gene_info$symbol[match(res2$variable, res$gene_info$gene)]
          left_pct <- (cpm_hover$x - cpm_hover$domain$left)/(cpm_hover$domain$right - cpm_hover$domain$left)
          top_pct <- (cpm_hover$domain$top - cpm_hover$y)/(cpm_hover$domain$top - cpm_hover$domain$bottom)
          left_px <- cpm_hover$range$left + left_pct * (cpm_hover$range$right - cpm_hover$range$left)
          top_px <- cpm_hover$range$top + top_pct * (cpm_hover$range$bottom - cpm_hover$range$top)
          style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                          "left:", left_px + 2, "px; top:", top_px + 2, "px;")
          wellPanel(
            style = style,
            p(HTML(paste0("<b> Sample: </b>", res2$sample, "<br/>",
                          "<b> ID: </b>", res2$variable, "<br/>",
                          "<b> Gene: </b>", res2$genename, "<br/>",
                          "<b> CPM: </b>", round(res2$value, 3), "<br/>")))
          )
        }
      }
    })
    
  }
  
  shinyApp(ui = p_layout, server = server_function, enableBookmarking = "server")
}

print(summary_app(res))
