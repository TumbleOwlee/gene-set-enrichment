#!/usr/bin/env Rscript

#' 
#' Maintainer: David Loewe <49597367+TumbleOwlee@users.noreply.github.com>
#' License: MIT
#' 

# Default Configuration
config <- data.frame(
    # Enable if debugging
    quiet = FALSE,
    # Input annotated file
    annotated = './private/out.emapper.annotations'
)

#' Generate the geneid to go/kegg mappings by expanding the table based on
#' elements in 'GOs' and 'KEGG_ko'.
main <- function() {
    # Initialize and install dependencies.
    ifh.info('Check and install missing dependencies...')
    ifh.init(c('numbers'), quiet = TRUE)

    # Setup CLI argument parser
    option_list = list(
        make_option(c('--annotated'), type = 'character', default = config$annotated,
                    help = 'eggNog-mapper annotation output file (e.g. out.emapper.annotations).', metavar = 'FILE'),
        make_option(c('--quiet'), action = 'store_true', default = config$quiet, help = 'Suppress all output.')
    )
    opt_parser = OptionParser(description = 'Generates GENEID-GO and GENEID-KEGG_ko mappings based on eggNog-mapper output.', option_list = option_list)
    opt = ifh.parse_args(opt_parser)

    # Read eggNOG-mapper annotation table
    ifh.info('Read all inputs...')
    annotated.df <- ifh.table.import(opt$annotated, header = FALSE, sep = '\t', comment.char = '#', strip.white = TRUE, quote = '')
    colnames(annotated.df) <- ifh.file.get.line(path = opt$annotated, prefix = '#', comment.prefix = '##', sep = '\t')
    
    # Rename 'query' to 'GENEID' for GSE
    colnames(annotated.df)[colnames(annotated.df) == 'query'] <- 'GENEID'
    ifh.success(' Annotated inputs loaded.')

    # Expand the table based on 'GOs' and export to CSV
    geneid.to.go <- ifh.table.expand(annotated.df, by = 'GOs', sep = '\t', split = ',', limit.to = c('GENEID'))
    ifh.table.export(geneid.to.go, 'geneid-to-go.csv')
    
    # Expand the table based on 'KEGG_ko' and export to CSV
    geneid.to.kegg <- ifh.table.expand(annotated.df, by = 'KEGG_ko', sep = '\t', split = ',', limit.to = c('GENEID'), remove.prefix = 'ko:', na = '-')
    ifh.table.export(geneid.to.kegg, 'geneid-to-kegg.csv')
}

# Load library 'ifhutil'
suppressMessages({library(ifhutil)})

# Execute `main()` with prettified log
ifh.run({main()})

###############################################################################
