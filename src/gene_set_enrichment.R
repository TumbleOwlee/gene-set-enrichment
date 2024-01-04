#!/usr/bin/env Rscript

#'
#' Maintainer: David Loewe <49597367+TumbleOwlee@users.noreply.github.com>
#' License: MIT
#'

# Setup default configuration
config <- data.frame(
    # Prevent any output
    quiet = TRUE,
    # Organism annotation database
    organism.db = "org.Dm.eg.db",
    # Input file
    input.file = "./input/drosphila_example_de.csv",
    # Keytype: source of annotation
    key.type = "ENSEMBL",
    # ont: "BP", "MF", "CC", "ALL"
    ont = "CC",
    # Number of permutations
    num.perm = 10000,
    # Minimum number of genes in set
    min.gs.size = 3,
    # Maximum number of genes in set
    max.gs.size = 800,
    # pvalue cutoff
    pvalue.cutoff = 0.05,
    # padjust method: "holm", "hochberg", "hommel", "bonferroni",
    #                 "BH", "BY", "fdr", "none"
    padjust.method = "none",
    # KEGG organism code, see https://www.genome.jp/kegg/catalog/org_list.html
    kegg.code = "dme",
    # KEGG key type
    kegg.key.type = "ncbi-geneid",
    # GO Map
    go.map = "",
    # KEGG Map
    kegg.map = "",
    # Run Gene Ontology
    run.go = TRUE,
    # Run KEGG
    run.kegg = TRUE,
    # Count of categories to plot
    show.categories = 10000,
    # Pathway ID
    pathway.id = "ecok00010",
    pathway.gene.idtype = "KEGG"
)

# Change configuration for e.coli K12 based on eggNOG-mapper output
if (TRUE) {
    config$organism.db <- "org.EcK12.eg.db"
    config$go.map <- "geneid-to-go.csv"
    config$go.toType <- "GENENAME"
    config$ont <- "BP"
    config$kegg.map <- "geneid-to-kegg.csv"
    config$input.file <- "./private/L12VLB12_MGaccession.csv"
    config$key.type <- "GO"
    config$kegg.key.type <- "kegg"
    config$kegg.code <- "ecok"
    config$run.kegg <- TRUE
    config$run.go <- TRUE
    config$pathway.id = "ecok00010",
    config$pathway.gene.idtype = "KEGG"
}

#' Performs GO and KEGG gene set enrichment analysis
main <- function() {
    # Initialize and install dependencies
    ifh.info("Check and install missing dependencies...")
    ifh.init(c("clusterProfiler", "pathview", "enrichplot", "gridExtra",
               "DOSE", "ggridges", "ggplot2", "tools", "annotate", "GO.db"), quiet = config$quiet)

    # Create CLI argument parser
    option_list = list(
        make_option(c("--file"),
                    type = "character",
                    default = config$input.file,
                    help = "Input file for gene enrichment analysis.",
                    metavar = "FILE"),
        make_option(c("--organism_annotation"),
                    type = "character",
                    default = config$organism.db,
                    help = "Input organism to use for gene enrichment analysis.",
                    metavar = "ANNOTATION"),
        make_option(c("--no_gene_ontology"),
                    action = "store_false",
                    default = !config$run.go,
                    help = "Disbale GO run."),
        make_option(c("--no_kyoto_genes"),
                    action = "store_false",
                    default = !config$run.kegg,
                    help = "Disbale KEGG run."),
        make_option(c("--go_map"),
                    type = "character",
                    default = config$go.map,
                    help = "Specify custom GENEID-to-GO mapping file.",
                    metavar = "FILE"),
        make_option(c("--go_toType"),
                    type = "character",
                    default = config$go.toType,
                    help = "Specify GO to type mapping.",
                    metavar = "TYPE"),
        make_option(c("--kegg_map"),
                    type = "character",
                    default = config$kegg.map,
                    help = "Specify custom GENEID-to-KEGG mapping file.",
                    metavar = "FILE"),
        make_option(c("--quiet"),
                    action = "store_true",
                    default = config$quiet,
                    help = "Suppress most of the logging."),
        make_option(c("--list_keytypes"),
                    action = "store_true",
                    default = FALSE,
                    help = "List options for source of annotation."),
        make_option(c("--keytype"),
                    type = "character",
                    default = config$key.type,
                    help = "Source of annotation to use for processing.",
                    metavar = "KEYTYPE"),
        make_option(c("--ont"),
                    type = "character",
                    default = config$ont,
                    help = "The ont value to use.",
                    metavar = "ONT"),
        make_option(c("--min_gs_size"),
                    type = "integer",
                    default = config$min.gs.size,
                    help = "Minimum number of genes in set.",
                    metavar = "MIN"),
        make_option(c("--max_gs_size"),
                    type = "integer",
                    default = config$max.gs.size,
                    help = "Maximum number of genes in set.",
                    metavar = "MAX"),
        make_option(c("--pvalue_cutoff"),
                    type = "double",
                    default = config$pvalue.cutoff,
                    help = "The pvalue cutoff.",
                    metavar = "CUTOFF"),
        make_option(c("--padjust_method"),
                    type = "character",
                    default = config$padjust.method,
                    help = "The padjust method name.",
                    metavar = "METHOD"),
        make_option(c("--kegg_code"),
                    type = "character",
                    default = config$kegg.code,
                    help = "The kegg organism code.",
                    metavar = "KEGG"),
        make_option(c("--kegg_keytype"),
                    type = "character",
                    default = config$kegg.key.type,
                    help = "The kegg key type.",
                    metavar = "KEY"),
        make_option(c("--show_categories"),
                    type = "integer",
                    default = config$show.categories,
                    help = "The count of categories to plot.",
                    metavar = "N"),
        make_option(c("--kegg_pathway_id"),
                    type = "character",
                    default = config$pathway.id,
                    help = "The kegg pathway id.",
                    metavar = "ID"),
        make_option(c("--kegg_pathway_gene_idtype"),
                    type = "character",
                    default = config$pathway.gene.idtype,
                    help = "The kegg pathway gene id type.",
                    metavar = "TYPE")
    )
    opt_parser = OptionParser(option_list = option_list)
    opt = ifh.parse_args(opt_parser)

    # Install/Load organism annotations
    ifh.info("Load organism annotation...")
    ifh.install(pkg = opt$organism_annotation, quiet = opt$quiet)
    ifh.success("Organism annotation loaded.")

    # Print supported key types of the organism annotations
    if (opt$list_keytypes) {
        ifh.info("Possible keytypes for", opt$organism_annotation, ":")
        ifh.info(keytypes(get(opt$organism_annotation)))
        return()
    }

    # Create cache directory if not present
    ifh.dir.create_if("./.cache")

    # Create GO cache file name
    go.geneid.map.file = ifh.create.path("./.cache/",
                                         ifh.create.filename(opt$organism_annotation, file_path_sans_ext(basename(opt$go_map))),
                                         "_go.RData")
    # Create KEGG cache file name
    kegg.geneid.map.file = ifh.create.path("./.cache/",
                                           ifh.create.filename(opt$organism_annotation,
                                                               file_path_sans_ext(basename(opt$kegg_map))),
                                           "_kegg.RData")
    # Create GO cache file path
    go.cache.file = ifh.create.path("./.cache/",
                                    ifh.create.filename(opt$organism_annotation, file_path_sans_ext(basename(opt$file))),
                                    "_go.RData")
    # Create KEGG cache file path
    kegg.cache.file = ifh.create.path("./.cache/",
                                      ifh.create.filename(opt$organism_annotation, file_path_sans_ext(basename(opt$file))),
                                      "_kegg.RData")

    # By default no calculation are executed if cache files exist
    run_calculations = FALSE

    # Load cache file or import go-geneid mapping frame
    go.geneid.df <- NA
    if (opt$go_map!= "") {
        if (ifh.cache.loadable(go.geneid.map.file)) {
            ifh.cache.load(go.geneid.map.file)
        } else {
            ifh.info("Import GO map from", opt$go_map)
            go.geneid.df <- ifh.table.import(opt$go_map, header = TRUE, sep = "\t", comment.char = "#")
            ifh.step("Create cache file.")
            ifh.cache.save(go.geneid.df, file = go.geneid.map.file)
            ifh.success("GO map imported.")
        }
    }

    # Load cache file or import KEGG-geneid mapping frame
    kegg.geneid.df <- NA
    if (opt$kegg_map != "") {
        if (ifh.cache.loadable(kegg.geneid.map.file)) {
            ifh.cache.load(kegg.geneid.map.file)
        } else {
            ifh.info("Import KEGG map from", opt$kegg_map)
            kegg.geneid.df <- ifh.table.import(opt$kegg_map, header = TRUE, sep = "\t", comment.char = "#")
            ifh.step("Create cache file.")
            ifh.cache.save(kegg.geneid.df, file = kegg.geneid.map.file)
            ifh.success("KEGG map imported.")
        }
    }

    # Load GO GSE cache file if present and requested
    if (!opt$no_gene_ontology) {
        if (ifh.cache.loadable(go.cache.file)) {
            ifh.cache.load(go.cache.file)
        } else {
            run_calculations = TRUE
        }
    }
    # Load KEGG GSE cache file if present and requested
    if (!opt$no_kyoto_genes) {
        if (ifh.cache.loadable(kegg.cache.file)) {
            ifh.cache.load(kegg.cache.file)
        } else {
            run_calculations = TRUE
        }
    }

    # If no cached data is present or user requested recalculation
    if (run_calculations) {
        # Read DESeq2 data
        ifh.step("Load CSV data from", opt$file, "...")
        deseq2.df = ifh.table.import(opt$file, header = TRUE, sep = ",")
        ifh.success("Data loaded.")

        ifh.step("Prepare gene list...")
        orig.gene.list <- deseq2.df$log2FoldChange
        names(orig.gene.list) <- deseq2.df$X
        gene.list <- na.omit(orig.gene.list)
        gene.list = sort(gene.list, decreasing = TRUE)
        ifh.success("Gene list prepared.")

        # Run GO GSE with organism database or custom GO-GENEID mapping
        if (!opt$no_gene_ontology) {
            ifh.step("Calculate GO gene set enrichment...")
            ifh.quiet(opt$quiet)
            go.gse = NA
            go.gene.list = NA

            if (opt$go_map != "") {
                go.supported.ids <- select(org.EcK12.eg.db, keys(org.EcK12.eg.db, "GO"), "ENTREZID", "GO")[,'GO']
                go.ids.in.ont <- filterGOByOntology(go.supported.ids, ontology = c(opt$ont))
                go.supported.ids <- go.supported.ids[go.ids.in.ont == TRUE]
                go.term.to.gene <- go.geneid.df[go.geneid.df$GO %in% go.supported.ids,]
                go.term.to.name <- go2term(go.geneid.df[,1])
                names(go.term.to.name) <- c('GO', 'NAME')

                go.gene.list <- gene.list

                go.gse = suppressMessages({
                    GSEA(geneList = go.gene.list,
                         exponent = 1,
                         minGSSize = opt$min_gs_size,
                         maxGSSize = opt$max_gs_size,
                         pvalueCutoff = opt$pvalue_cutoff,
                         pAdjustMethod = opt$padjust_method,
                         verbose = FALSE,
                         eps = 0,
                         TERM2GENE = go.term.to.gene,
                         TERM2NAME = go.term.to.name,
                         by = "fgsea"
                         )
                })
                
                convert_go_to_pref_name <- function(val, conversion_table) {
                    values <- conversion_table$PreferredName[conversion_table$GENEID %in% val]
                    values <- unique(values)
                    return(paste(values, collapse = ',', sep = ','))
                }
                
                #mapping
                convert_list_go_to_pref_name <- function(list, conversion_table) {
                    vec <- unlist(strsplit(list, "/"))
                    gene_names <- sapply(vec, convert_go_to_pref_name, conversion_table = conversion_table)
                    return(paste(gene_names, collapse = "/"))
                }
                
                #extract the dataframe from enrichKEGG object
                temp_df = go.gse@result
                
                #replace the IDs with sapply
                temp_df$core_enrichment = sapply(temp_df$core_enrichment, convert_list_go_to_pref_name, conversion_table = go.geneid.df)
                
                #modify the result table
                go.gse@result = temp_df
                
                for (key in names(go.gse@geneSets)) {
                    go.gse@geneSets[key] <- sapply(go.gse@geneSets[key], convert_go_to_pref_name, conversion_table = go.geneid.df)
                }
                
                names(go.gse@geneList) <- sapply(names(go.gse@geneList), convert_go_to_pref_name, conversion_table = go.geneid.df)
                names(go.gene.list) <- sapply(names(go.gene.list), convert_go_to_pref_name, conversion_table = go.geneid.df)
            } else {
                go.gene.list <- gene.list
                go.gse = suppressMessages({
                    gseGO(geneList = go.gene.list,
                          ont = opt$ont,
                          keyType = opt$keytype,
                          #nPerm = 10000,
                          minGSSize = opt$min_gs_size,
                          maxGSSize = opt$max_gs_size,
                          pvalueCutoff = opt$pvalue_cutoff,
                          verbose = FALSE,
                          OrgDb = opt$organism_annotation,
                          pAdjustMethod = opt$padjust_method)})
            }
            ifh.quiet(FALSE)
            ifh.step("Calculations for GO GSE finished.")

            # Store cached results
            ifh.info("Save results in cache file", go.cache.file)
            ifh.cache.save(go.gse, go.gene.list, file = go.cache.file)
        }

        # Run KEGG GSE with organism database or custom KEGG_ko-GENEID mapping
        if (!opt$no_kyoto_genes) {
            ifh.step("Prepare KEGG gene list...")

            ifh.quiet(opt$quiet)
            if (opt$kegg_map != "") {
                kegg.df.unique.ids = kegg.geneid.df[!duplicated(kegg.geneid.df[c("GENEID")]),]
                deseq2.df = deseq2.df[deseq2.df$X %in% kegg.df.unique.ids[["GENEID"]],]
                kegg.df.unique.ids = kegg.df.unique.ids[kegg.df.unique.ids$GENEID %in% deseq2.df$X,]
                deseq2.df$Y = kegg.df.unique.ids$KEGG_ko
                
                kegg.gene.list <- deseq2.df$log2FoldChange
                names(kegg.gene.list) <- deseq2.df$Y
                kegg.gene.list <- na.omit(kegg.gene.list)
                kegg.gene.list = sort(kegg.gene.list, decreasing = TRUE)
                ifh.success("KEGG gene list prepared.")

                ifh.step("Calculate KEGG gene set enrichment...")
                kegg.gse <- gseKEGG(geneList = kegg.gene.list,
                                    organism = "ko",
                                    #nPerm = 10000,
                                    minGSSize = opt$min_gs_size,
                                    maxGSSize = opt$max_gs_size,
                                    pvalueCutoff = opt$pvalue_cutoff,
                                    pAdjustMethod = opt$padjust_method,
                                    keyType = opt$kegg_keytype,
                                    eps = 0)
                
                convert_kegg_to_pref_name <- function(val, conversion_table) {
                    values <- conversion_table$PreferredName[conversion_table$KEGG_ko %in% val]
                    values <- unique(values) 
                    return(paste(values, collapse = ',', sep = ','))
                }
                
                #mapping
                convert_list_kegg_to_pref_name <- function(list, conversion_table) {
                    vec <- unlist(strsplit(list, "/"))
                    gene_names <- sapply(vec, convert_kegg_to_pref_name, conversion_table = conversion_table)
                    return(paste(gene_names, collapse = "/"))
                }
                
                #extract the dataframe from enrichKEGG object
                temp_df = kegg.gse@result
                
                #replace the IDs with sapply
                temp_df$core_enrichment = sapply(temp_df$core_enrichment, convert_list_kegg_to_pref_name, conversion_table = kegg.geneid.df)
                
                #modify the result table
                kegg.gse@result = temp_df
                
                for (key in names(kegg.gse@geneSets)) {
                    kegg.gse@geneSets[key] <- sapply(kegg.gse@geneSets[key], convert_kegg_to_pref_name, conversion_table = kegg.geneid.df)
                }
                
                names(kegg.gse@geneList) <- sapply(names(kegg.gse@geneList), convert_kegg_to_pref_name, conversion_table = kegg.geneid.df)
            } else {
                ids <- bitr(names(orig.gene.list), fromType = opt$keytype, toType = "ENTREZID", OrgDb = opt$organism_annotation)

                kegg.df.unique.ids = ids[!duplicated(ids[c(opt$keytype)]),]
                deseq2.df = deseq2.df[deseq2.df$X %in% kegg.df.unique.ids[[opt$keytype]],]
                deseq2.df$Y = kegg.df.unique.ids$ENTREZID

                kegg.gene.list <- deseq2.df$log2FoldChange
                names(kegg.gene.list) <- deseq2.df$Y
                kegg.gene.list <- na.omit(kegg.gene.list)
                kegg.gene.list = sort(kegg.gene.list, decreasing = TRUE)
                ifh.success("KEGG gene list prepared.")

                ifh.step("Calculate KEGG gene set enrichment...")
                kegg.gse <- gseKEGG(geneList = kegg.gene.list,
                                    organism = opt$kegg_code,
                                    #nPerm = 10000,
                                    minGSSize = opt$min_gs_size,
                                    maxGSSize = opt$max_gs_size,
                                    pvalueCutoff = opt$pvalue_cutoff,
                                    pAdjustMethod = opt$padjust_method,
                                    keyType = opt$kegg_keytype)
            }
            ifh.quiet(FALSE)
            ifh.step("Calculations for KEGG GSE finished.")

            # Store cached results
            ifh.info("Save results in cache file", kegg.cache.file)
            ifh.cache.save(kegg.gse, kegg.gene.list, gene.list, file = kegg.cache.file)
        }
    }

    # Plot GO GSE results
    if (!opt$no_gene_ontology) {
        ifh.step("Create GO dotplot...")
        p <- dotplot(go.gse, showCategory = opt$show_categories, split = ".sign", font.size=8) + facet_grid(.~.sign)
        print(p)
        ifh.success("Dotplot created.")

        # ifh.step("Create GO emapplot...")
        # go.gse2 <- pairwise_termsim(go.gse, method="JC")
        # p <- emapplot(go.gse2, showCategory = 10, cex.params = list(category_label = 0.5))
        # print(p)
        # ifh.success("Emapplot created.")
        #
        ifh.step("Create GO cnetplot...")
        p <- cnetplot(go.gse, categorySize="pvalue", color.params = list(foldChange=go.gene.list), showCategory = opt$show_categories, colorEdge = TRUE, cex.params = list(category_label = 0.7, gene_label = 0.5))
        print(p)
        ifh.success("Cnetplot created.")
        
        ifh.step("Create GO heatplot...")
        p <- heatplot(go.gse, foldChange=go.gene.list, showCategory = opt$show_categories)
        print(p)
        ifh.success("Heatplot created.")
        
        #
        # ifh.step("Create GO ridgeplot...")
        # p <- ridgeplot(go.gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size=5), axis.text.x = element_text(size=8))
        # print(p)
        # ifh.success("Ridgeplot created.")
        #
        # ifh.step("Create GO gseaplot...")
        # p <- gseaplot(go.gse, by = "all", title = go.gse$Description[1], geneSetID = 1) + theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=8))
        # print(p)
        # ifh.success("Gseaplot created.")
        #
        # ifh.step("Create GO pmcplot...")
        # terms <- go.gse$Description[1:3]
        # p <- pmcplot(terms, 2010:2023, proportion = FALSE) + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        # print(p)
        # ifh.success("Pmcplot created.")
    }

    # Plot KEGG GSE results
    if (!opt$no_kyoto_genes) {
        ifh.step("Create KEGG dotplot...")
        p <- dotplot(kegg.gse, showCategory = opt$show_categories, title = "Enriched Pathways", split = ".sign", font.size=8) + facet_grid(.~.sign)
        print(p)
        ifh.success("Dotplot created.")

        # ifh.step("Create KEGG emapplot...")
        # kegg.gse2 <- pairwise_termsim(kegg.gse, method="JC")
        # p <- emapplot(kegg.gse2, cex.params = list(category_label = 0.5))
        # print(p)
        # ifh.success("Emapplot created.")
        #
        ifh.step("Create KEGG cnetplot...")
        
        p <- cnetplot(kegg.gse, categorySize="pvalue", color.params = list(foldChange=kegg.gse@geneList), showCategory = opt$show_categories, colorEdge = TRUE, cex.params = list(category_label = 0.7, gene_label = 0.5))
        print(p)
        ifh.success("Cnetplot created.")
        
        ifh.step("Create KEGG heatplot...")
        p <- heatplot(kegg.gse, foldChange=kegg.gse@geneList, showCategory = opt$show_categories)
        print(p)
        ifh.success("Heatplot created.")
        
        #
        # ifh.step("Create KEGG ridgeplot...")
        # p <- ridgeplot(kegg.gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size=5), axis.text.x = element_text(size=8))
        # print(p)
        # ifh.success("Ridgeplot created.")
        #
        # ifh.step("Create KEGG gseaplot...")
        # p <- gseaplot(kegg.gse, by = "all", title = kegg.gse$Description[1], geneSetID = 1) + theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=8))
        # print(p)
        # ifh.success("Gseaplot created.")


        # Display KEGG pathview (currently disabled)
    
        show.pathview <- function(..., save = FALSE)
        {
            ifh.info("Create image..")
            msg <- capture.output(pathview::pathview(...), type = "message")
            msg <- grep("image file", msg, value = T)
            filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])

              img <- png::readPNG(filename)
              grid::grid.raster(img)

              filename2 = gsub("\\.pathview", "", filename)
              img2 <- png::readPNG(filename2)
              grid::grid.raster(img2)

              filename3 = gsub("\\.png", ".xml", filename2)

              if(!save) {
                  invisible(file.remove(filename))
                  invisible(file.remove(filename2))
                  invisible(file.remove(filename3))
              }
        }
        
        ifh.step("Create KEGG pathview...")
        pathview(cpd.data=kegg.gene.list, pathway.id = opt$kegg_pathway_id, species = opt$kegg_code, gene.idtype = opt$kegg_pathway_gene_idtype)
        ifh.success("Pathview created.")
    }
}

################################### EXECUTE ###################################

suppressMessages({library(ifhutil)})

#ifh.run({main()})
main()

###############################################################################
