pkgname <- "scicomptools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scicomptools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("drive_toc")
### * drive_toc

flush(stderr()); flush(stdout())

### Name: drive_toc
### Title: Identify all Folders within Specified Google Drive Folder
### Aliases: drive_toc

### ** Examples

 
## Not run: 
##D # Supply a single Google Drive folder link to identify all its sub-folders 
##D drive_toc(url = googledrive::as_id("https://drive.google.com/drive/u/0/folders/your-folder"))
## End(Not run)




cleanEx()
nameEx("read_xl_format")
### * read_xl_format

flush(stderr()); flush(stdout())

### Name: read_xl_format
### Title: Read Formatting of All Sheets in an Excel Workbook
### Aliases: read_xl_format

### ** Examples

# Identify the formatting of every cell in all sheets of an Excel file
read_xl_format(file_name = system.file("extdata", "excel_book.xlsx", package = "scicomptools"))




cleanEx()
nameEx("read_xl_sheets")
### * read_xl_sheets

flush(stderr()); flush(stdout())

### Name: read_xl_sheets
### Title: Read All Sheets from an Excel Workbook
### Aliases: read_xl_sheets

### ** Examples

# Read in each sheet as an element in a list
read_xl_sheets(file_name = system.file("extdata", "excel_book.xlsx", package = "scicomptools"))




cleanEx()
nameEx("stat_extract")
### * stat_extract

flush(stderr()); flush(stdout())

### Name: stat_extract
### Title: Extract Summary Statistics from Model Fit Object
### Aliases: stat_extract

### ** Examples

# Create some example data
x <- c(3.5, 2.1, 7.5, 5.6, 3.3, 6.0, 5.6)
y <- c(2.3, 4.7, 7.8, 9.1, 4.5, 3.6, 5.1)

# Fit a linear model
mod <- lm(y ~ x)

# Extract the relevant information
stat_extract(mod_fit = mod)




cleanEx()
nameEx("token_check")
### * token_check

flush(stderr()); flush(stdout())

### Name: token_check
### Title: Check Token Status
### Aliases: token_check

### ** Examples

## Not run: 
##D # Check whether a GitHub token is attached or not
##D token_check(api = "github", secret = TRUE)
## End(Not run)
## Not run: 
##D # Check whether a Qualtrics token is attached or not
##D token_check(api = "qualtrics", secret = TRUE)
## End(Not run)



cleanEx()
nameEx("wd_loc")
### * wd_loc

flush(stderr()); flush(stdout())

### Name: wd_loc
### Title: Define Local or Remote Working Directories
### Aliases: wd_loc

### ** Examples

# Set two working directory paths to toggle between

# If you are working in your local computer, set `local` to "TRUE"
wd_loc(local = TRUE,
       local_path = file.path("local path"),
       remote_path = file.path("path on server"))
       
# If you are working in a remote server, set `local` to "FALSE"
wd_loc(local = FALSE,
       local_path = file.path("local path"),
       remote_path = file.path("path on server"))
      



cleanEx()
nameEx("word_cloud_plot")
### * word_cloud_plot

flush(stderr()); flush(stdout())

### Name: word_cloud_plot
### Title: Text Mine a Given Column and Create a Word Cloud
### Aliases: word_cloud_plot

### ** Examples

# Create a dataframe containing some example text
text <- data.frame(article_num = 1:6,
                   article_title = c("Why pigeons are the best birds",
                                     "10 ways to show your pet budgie love",
                                     "Should you feed ducks at the park?",
                                     "Locations and tips for birdwatching",
                                     "How to tell which pet bird is right for you",
                                     "Do birds make good pets?"))
                                     
# Prepare the dataframe for word cloud plotting              
word_cloud_prep(data = text, text_column = "article_title")

# Plot the word cloud
word_cloud_plot(data = text, text_column = "article_title")




cleanEx()
nameEx("word_cloud_prep")
### * word_cloud_prep

flush(stderr()); flush(stdout())

### Name: word_cloud_prep
### Title: Perform Text Mining of a Given Column
### Aliases: word_cloud_prep

### ** Examples

# Create a dataframe containing some example text
text <- data.frame(article_num = 1:6,
                   article_title = c("Why pigeons are the best birds",
                                     "10 ways to show your pet budgie love",
                                     "Should you feed ducks at the park?",
                                     "Locations and tips for birdwatching",
                                     "How to tell which pet bird is right for you",
                                     "Do birds make good pets?"))
                                     
# Prepare the dataframe for word cloud plotting              
word_cloud_prep(data = text, text_column = "article_title")

# Plot the word cloud
word_cloud_plot(data = text, text_column = "article_title")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
