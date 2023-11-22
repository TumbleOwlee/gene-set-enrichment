#'
#' Maintainer: David Loewe <49597367+TumbleOwlee@users.noreply.github.com>
#' License: MIT
#'

#' Print a message containing the concatenated string of all arguments.
#'
#' @description
#' This function prints the given message with the given `prefix` if provided.
#' It is mainly used by `ifh.info()`, `ifh.step()`, `ifh.success()` and `ifh.error()`.
#'
#' @param ...    arguments passed down to paste(..)
#' @param prefix A prefix to add to the message
#' @examples
#' ifh.message("The number is", 10)
#' ifh.message("Add prefix to mark the message.", prefix = "[!]")
#'
#' @export
ifh.message <- function(..., prefix = "") {
    if (prefix == "") {
        message(paste(..., sep = " ", collapse = " "))
    } else {
        message(paste(prefix, ..., sep = " ", collapse = " "))
    }
}

#' Prints the arguments as concatenated string and marked as info.
#'
#' @description
#' This function prints the given message and adds the log level prefix
#' `[i] INFO:`. After the call to `ifh.init()`, this function is replaced by
#' `cli_alert_info()` by the 'cli' package for prettier messages.
#'
#' @param ...   arguments passed down to paste(..)
#' @examples
#' ifh.info("Print some info with number", 10)
#'
#' @export
ifh.info <- function(...) {
    ifh.message(paste(..., sep = " ", collapse = " "), prefix = "[i] INFO:")
}

#' Prints the arguments as concatenated string and marked as warning.
#'
#' @description
#' This function prints the given message and adds the log level prefix
#' `[!] WARN:`. After the call to `ifh.init()`, this function is replaced by
#' `cli_alert_warning()` by the 'cli' package for prettier messages.
#'
#' @param ...   arguments passed down to paste(..)
#' @examples
#' ifh.warn("Print some step with number", 10)
#'
#' @export
ifh.warn <- function(...) {
    ifh.message(paste(..., sep = " ", collapse = " "), prefix = "[!] WARN:")
}

#' Prints the arguments as concatenated string and marked as step.
#'
#' @description
#' This function prints the given message and adds the log level prefix
#' `[-] STEP:`. After the call to `ifh.init()`, this function is replaced by
#' `cli_alert()` by the 'cli' package for prettier messages.
#'
#' @param ...   arguments passed down to paste(..)
#' @examples
#' ifh.step("Print some step with number", 10)
#'
#' @export
ifh.step <- function(...) {
    ifh.message(paste(..., sep = " ", collapse = " "), prefix = "[-] STEP:")
}

#' Prints the arguments as concatenated string and marked as success.
#'
#' @description
#' This function prints the given message and adds the log level prefix
#' `[!] SUCCESS:`. After the call to `ifh.init()`, this function is replaced by
#' `cli_alert_success()` by the 'cli' package for prettier messages.
#'
#' @param ...   arguments passed down to paste(..)
#' @examples
#' ifh.success("Print a success message with number", 10)
#'
#' @export
ifh.success <- function(...) {
    ifh.message(paste(..., sep = " ", collapse = " "), prefix = "[!] SUCCESS:")
}

#' Prints the arguments as concatenated string and marked as error
#'
#' @description
#' This function prints the given message and adds the log level prefix
#' `[e] ERROR:`. After the call to `ifh.init()`, this function is replaced by
#' `cli_alert_danger()` by the 'cli' package for prettier messages.
#'
#' @param ...   arguments passed down to paste(..)
#' @examples
#' ifh.error("Print an error message with number", 10)
#'
#' @export
ifh.error <- function(...) {
    ifh.message(paste(..., sep = " ", collapse = " "), prefix = "[e] ERROR:")
}

#' Execute the given expression either quietly or not
#'
#' @description
#' This function is used to execute any expression and providing the option
#' to suppress all output by utilizing a flag. It is mainly used by CLI
#' applications that provide an `--quiet` flag to switch output.
#'
#' @param expr  expression to execute
#' @param quiet if true, all output of expression is suppressed
#' @examples
#' ifh.exec({install.packages('cli')}, quiet = TRUE)
#'
#' @export
ifh.exec <- function(expr, quiet = FALSE) {
    if (quiet == TRUE) {
        suppressMessages(expr)
    } else {
        eval(expr)
    }
}

#' Install and load the given package using `BiocManager`, usable after `ifh.init(..)`
#'
#' @description
#' This function allows the installation of a package. It is only callable
#' without an error after `ifh.init()` is called. It can be used to install and load
#' additional packages after the `ifh.init()` was called.
#'
#' @param pkg   package name to install
#' @param ...   arguments passed down to install(..)
#' @param quiet suppress all output or not
#' @examples
#' ifh.install('cli', quiet = FALSE)
#'
ifh.install <- function(pkg, ..., quiet = FALSE) {
    stop("Not initialized. Make sure 'ifh.init(..)' was called!")
}

#' Initializes the project by setting up logging and libraries
#'
#' @description
#' This function installs 'cli', 'progress' and 'optparse' by default
#' for pretty output. Furthermore any list of package names provided as
#' `packages` are installed and loaded. Also the flag is provided to
#' suppress all output of the installation and load calls.
#'
#' @param packages  vector of all packages to install
#' @param quiet     suppress all output or not
#'
#' @examples
#' ifh.init(c('somePackage', 'otherPackage'), quiet = TRUE)
#'
#' @export
ifh.init <- function(packages, quiet = FALSE) {
    # Install BiocManager is not present
    if (!require("BiocManager", quietly = quiet))
        install.packages("BiocManager", repos = "http://cran.rstudio.com/", dependencies = TRUE, quiet = quiet)

    # Get list of installed packages to skip install calls for
    # already present packages.
    installed <- rownames(installed.packages())

    # Initialize the install function to use BiocManager and install
    # packages only not present already.
    .GlobalEnv$ifh.install <- function(pkg, ..., quiet = FALSE) {
        if (!(pkg %in% installed)) {
            ifh.exec({BiocManager::install(pkg = pkg, dependencies = TRUE, quiet = quiet, ...)}, quiet = quiet)
        }
        ifh.exec({library(pkg, character.only = TRUE, quietly = quiet)}, quiet = quiet)
    }

    # Install package 'cli' by default (skip if present in provided packages list)
    # This package is used for prettier log messages.
    if (!("cli" %in% packages)) {
        .GlobalEnv$ifh.install(pkg = "cli", quiet = quiet)
    }

    # Install package 'progress' by default (skip if present in provided packages list)
    # This package is used for displaying progress bars. This way the user gets feedback
    # and knows how much work is left.
    if (!("progress" %in% packages)) {
        .GlobalEnv$ifh.install(pkg = "progress", quiet = quiet)
    }

    # Install package 'optparse' by default (skip if present in provided packages list)
    # This package is used to provide a CLI to the user.
    if (!("optparse" %in% packages)) {
        .GlobalEnv$ifh.install(pkg = "optparse", quiet = quiet)
    }

    # Install all user requested packages.
    for (pkg in packages) {
        .GlobalEnv$ifh.install(pkg = pkg, quiet = quiet)
    }

    # Exchange the previous logging function with calls using package 'cli'
    .GlobalEnv$ifh.info <- function(...) { cli_alert_info(paste("", ..., sep = " ")) }
    .GlobalEnv$ifh.step <- function(...) { cli_alert(paste("", ..., sep = " ")) }
    .GlobalEnv$ifh.success <- function(...) { cli_alert_success(paste("", ..., sep = " ")) }
    .GlobalEnv$ifh.error <- function(...) { cli_alert_danger(paste("", ..., sep = " ")) }
    .GlobalEnv$ifh.warn <- function(...) { cli_alert_warning(paste("", ..., sep = " ")) }

    # Wrap optparse::parse_args() to print all set options. This allows the user
    # to see whether some requested option was recognized or not.
    .GlobalEnv$ifh.parse_args <- function(parser, ...) {
        opt <- optparse::parse_args(parser, ...)
        .GlobalEnv$ifh.info("Active options:")
        for (name in names(opt)) {
            .GlobalEnv$ifh.step(name, ":", opt[name])
        }
        return(opt)
    }
}

#' Executes the expression and adds pretty output of warnings and errors
#'
#' @description
#' This function wraps the given expression and adds handling of warning and
#' error messages. Especially it uses the 'cli' package to prettify the
#' warning and error messages.
#'
#' @param expr  expression to execute
#'
#' @export
ifh.run <- function(expr) {
    tryCatch(
        withCallingHandlers(
            expr,
            warning = function(e) {
                .GlobalEnv$ifh.warn(e$message)
            }
        ),
        error = function(e) {
            .GlobalEnv$ifh.error(e$message)
        }
    )
}

#' Retrieves user input in (non-)interactive mode
#'
#' @description
#' This function prompts the user for input. This wrapper is necessary
#' to provide the prompt functionality in non-interactive and interactive
#' mode.
#'
#' @param prompt    prompt to show to user
#'
#' @export
ifh.user.prompt <- function(prompt) {
    if (interactive()) {
        return(readline(prompt))
    } else {
        cat(prompt)
        return(readLines("stdin", n=1))
    }
}

#' Checks if file exists with given path
#'
#' @description
#' This function checks if the given file exists. It wraps the calls to
#' `file.exists()` and `dir.exists()` to only return TRUE if its a file and
#' not a directory.
#'
#' @param path  path to file to check
#' @return TRUE if file exists, FALSE otherwise
#'
#' @export
ifh.file.exists <- function(path) {
    return(file.exists(path) && !dir.exists(path))
}

#' Removes the given file
#'
#' @description
#' This function removes the given file if present. It returns TRUE
#' if it successfully removed the file. If the file doesn't exist,
#' it doesn't throw an error or warning but will return FALSE.
#'
#' @param path path to file for removal
#' @return TRUE if removed, FALSE if file doesn't exist
#'
#' @export
ifh.file.remove <- function(path) {
    if (ifh.file.exists(path)) {
        file.remove(path)
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Checks if directory exists
#'
#' @description
#' This function is currently a simple wrapper around `dir.exists()`.
#'
#' @param path  path to directory
#' @return TRUE if directory exists, FALSE otherwise
#'
#' @export
ifh.dir.exists <- function(path) {
    return(dir.exists(path))
}

#' Creates the given directory
#'
#' @description
#' This function is currently a simple wrapper around `dir.create()`.
#'
#' @param path  path to the directory
#'
#' @export
ifh.dir.create <- function(path) {
    dir.create(path)
}

#' Creates the directory if not already present
#'
#' @description
#' This function will check if the directory already exists. If so, it
#' returns `FALSE`. If not, the directory is created and `TRUE` is returned.
#'
#' @param path  path to the directory
#' @return TRUE if directory was created, FALSE if directory is already present
#'
#' @export
ifh.dir.create_if <- function(path) {
    if (ifh.dir.exists(path)) {
        return(FALSE)
    } else {
        ifh.dir.create(path)
        return(TRUE)
    }
}

#' Checks if the given cache file exists and user wants to load it
#'
#' @description
#' This function checks whether the cache file exists and if so, whether
#' the user wants to load it or not. If the cache file doesn't exist, FALSE
#' is returned. If the cache exists, the user is asked whether the cache
#' should be loaded or not. If the user decides to load it, TRUE is returned.
#' Else FALSE is returned.
#'
#' @param file   file path of the cache file
#' @param prompt Prompt to show to the user if cache exists
#' @return TRUE if cache exists and user wants to load it, FALSE otherwise
#'
#' @export
ifh.cache.loadable <- function(file, prompt = "Should the cache be loaded?") {
    if (ifh.file.exists(file)) {
        .GlobalEnv$ifh.info(paste(" Found cache file \"", file, "\"", sep = ""))
        input <- NA
        while (is.na(input) || (input != "y" && input != "Y" && input != "n" && input != "N")) {
            input <- ifh.user.prompt(paste("Cache exists.", prompt, "[Y]es [N]o ", sep = " "))
        }
        return(input == "y" || input == "Y")
    } else {
        return(FALSE)
    }
}

#' Save the given variables in a cache file
#'
#' @description
#' This functions saves all given variables in a cache file. Currently it
#' is a simple wrapper around `save()`.
#'
#' @param ...   arguments passed down to save(..)
#' @param file  output file name
#'
#' @export
ifh.cache.save <- function(..., file) {
    save(..., envir = parent.frame(), file = file)
}

#' Load the given cache file into environment
#'
#' @description
#' This function loads the given file into the environment. Currently it
#' is a simple wrapper around `load()`.
#'
#' @param file  file path to load
#'
#' @export
ifh.cache.load <- function(file) {
    load(file = file, envir = parent.frame())
}

#' Measures the performance of the given function call
#'
#' @description
#' This function measures the time the given function takes to be executed.
#' It prints the result via `ifh.info()`. This function can be used if you want
#' to determine the potential bottleneck that consumes a lot of CPU time.
#'
#' @param func  the function to measure
#' @param ...   arguments passed to func
#'
#' @return Return value of func
#'
#' @export
ifh.performance.measure <- function(func, ...) {
    start.time <- Sys.time()
    ret <- func(...)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    ifh.info(paste("Execution took ", time.taken))
    return(ret)
}

#' Creates n random strings
#'
#' @description
#' This function returns a list of n random generated string. It can be useful
#' in cases where you have to create temporary files that shouldn't clash with
#' existing files.
#'
#' @param n  number of random strings to create
#' @return vector of random strings
#'
#' @export
ifh.string.random <- function(n = 1) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#' Expands a table by duplicating a row for each element of a cell.
#'
#' @description
#' This function expands a table based on a given column containing a list of elements,
#' so that each row of the table only contains a single value in the column. It is done
#' by multiplying the row for each of the list elements. This is useful if you want to
#' create an 1:1 mapping from an 1:n mapping in each row.
#'
#' @param table          the data.frame to expand
#' @param by             the column name containing the list of elements
#' @param split          the delimiter of each element of the list
#' @param na             the NA value of the column
#' @param limit.to       column names to include in returned data.frame
#' @param remove.prefix  the prefix to remove from each list element before duplication
#' @return data.frame of all expanded rows with columns specified in 'by' and 'limit.to'
#'
#' @examples
#' ifh.table.expand(df, 'X', split = ',', na = '-', limit.to = c('Y'), remove.prefix = 'p:')
#'
#' @export
ifh.table.expand <- function(table, by, split = ",", na = "", limit.to = vector(mode="character"), remove.prefix = "") {
    sep = "\t"
    # Create progress bar to give user feedback
    # This is useful since the expansion can take some time based on the row count.
    pb <- progress_bar$new(total = nrow(table))
    pb$tick(0)

    .GlobalEnv$ifh.step(paste(nrow(table), "rows to expand."))

    # We create a cache file and print each expanded row into it.
    # This is useful since we can not predict the size of the result table
    # and adding a additional row to an existing data.frame multiple times
    # takes a lot of computational time since the data.frame has to reallocate
    # memory each time and copy the old memory into the new.
    # Loading a CSV file into a data.frame on the other hand is way faster.
    filename <- paste(".", ifh.string.random(1), ".csv", sep = "")
    output <- file(filename, open = "a")

    # Write header line
    writeLines(paste(by, paste(limit.to, collapse = sep, sep = sep), collapse = sep, sep = sep), output)

    for (i in 1:nrow(table)) {
        # Split list elements in given column
        values <- strsplit(table[i, by], split = split)[[1]]

        # For each element copy the row
        for (value in values) {
            # Remove the prefix if applicable
            if (value != na && remove.prefix != "") {
                value <- strsplit(value, split = remove.prefix)[[1]][2]
            }
            # Write new row
            if (is.na(value) || value == na) {
                writeLines(paste(NA, table[i, limit.to], sep = sep, collapse = sep), output)
            } else {
                writeLines(paste(value, table[i, limit.to], sep = sep, collapse = sep), output)
            }
        }

        # Inform user about progress
        pb$tick()
    }
    close(output)

    # Now read the expanded table
    data <- read.csv(filename, header = TRUE, sep = sep)
    # Remove the cache file
    file.remove(filename)
    # Inform user of expansion size
    .GlobalEnv$ifh.success(paste("All rows expanded. Total of ", nrow(data), "rows."))

    return(data)
}

#' Export the table as valid csv
#'
#' @description
#' This function exports the given table as valid CSV using the given separator.
#' It is used in cases where the `write.csv()` function won't work because it only allows
#' the comma separator.
#'
#' @param table     the data.frame to export
#' @param file      the output file name
#' @param row.names defines whether row.names are included
#' @param quote     defines whether quotes are included
#' @param sep       seperator to use as cell delimiter
#' @param ...       additional values to pass down to write.table(..)
#'
#' @export
ifh.table.export <- function(table, file, row.names = FALSE, quote = FALSE, sep = "\t", ...) {
    .GlobalEnv$ifh.info("Export table to", file)
    write.table(table, file = file, row.names = row.names, quote = quote, sep = sep, ...)
    .GlobalEnv$ifh.success("Export completed.")
}

#' Imports the given CSV file
#'
#' @description
#' This function is currently a simple wrapper around `read.csv()` with additional logging.
#'
#' @param file   the file path to import
#' @param sep    the separator of each cell
#' @param ...    additional parameters to pass down to read.csv(..)
#' @return the imported data.frame
#'
#' @export
ifh.table.import <- function(file, sep = "\t", ...) {
    .GlobalEnv$ifh.info("Import table from", file)
    table <- read.csv(file = file, sep = sep, ...)
    .GlobalEnv$ifh.success("Import completed.")
    return(table)
}

#' Retrieves the first line with the given prefix
#'
#' @description
#' This method can be used to retrieve the header of a table from a file
#' if the header has a special line prefix that prevents `read.csv(..)` from detecting
#' the header names.
#'
#' @param path            the path to the file
#' @param prefix          the prefix to detect
#' @param comment.prefix  the prefix of comments to ignore
#' @param sep             the separator of table cells
#' @return the first line matching the prefix (prefix removed) or none
#'
#' @export
ifh.file.get.line <- function(path, prefix, comment.prefix = "##", sep = "\t") {
    input <- file(path, "r")
    while(TRUE) {
        line <- readLines(input, 1)
        if (!startsWith(line, comment.prefix) && startsWith(line, prefix)) {
            split <- strsplit(substr(line, length(prefix) + 1, 100000000), "\t")[[1]]

            close(input)
            return(split)
        }
        if (length(line) == 0) break
    }
    return()
}

#' Creates a valid string path from the given arguments
#'
#' @description
#' This function is a simple wrapper around `paste(..., sep = '')`.
#'
#' @param ...   arguments passed down to paste(..)
#' @return valid path
#'
#' @export
ifh.create.path <- function(...) {
    return(paste(..., sep = ""))
}

#' Creates a valid file name form the given arguments
#'
#' @description
#' This function creates a valid file name from the given values by
#' concatenation using the separator '__' and replacement of '.' with '_'.
#'
#' @param ...   arguments passed down to paste(..)
#' @return valid file name
#'
#' @export
ifh.create.filename <- function(...) {
    combined <- paste(..., sep = "__")
    return(gsub("\\.", "_", combined))
}

#' Disables/Enables all output by redirecting to /dev/null
#'
#' @description
#' This function sets the sink to `/dev/null` or resets the sink
#' to the default based on the boolean value of quiet.
#'
#' @param quiet defines whether output is disabled or not
#'
#' @export
ifh.quiet <- function(quiet) {
    if (quiet == TRUE) {
        sink("/dev/null")
    } else {
        sink()
    }
}
