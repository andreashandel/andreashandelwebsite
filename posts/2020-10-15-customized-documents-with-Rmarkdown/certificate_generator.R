#####################################################
# automatically create customized certificates
#####################################################

# load needed packages
library('readr')
library('dplyr')
library('stringr')

# load data, read everything in as a string/character
df <- read.csv("student_data.csv",colClasses = 'character')

#set this to TRUE if you want to generate word output, FALSE for pdf
word_out = TRUE

# load either pdf or word certificate template
template <- ifelse(word_out, readr::read_file("certificate_template_word.Rmd"), readr::read_file("certificate_template_pdf.Rmd"))

#run through all students, generate personalized certificate for each
for (i in 1:nrow(df))
{

  #replace the placeholder words in the template with the student information
  current_cert <- template %>%
    str_replace("<<FIRSTNAME>>", df[i,'FirstName']) %>%
    str_replace("<<LASTNAME>>", df[i,'LastName']) %>%
    str_replace("<<SCORE>>", df[i,'Score'])

  #generate an output file name based on student name
  out_filename = paste(df[i,'LastName'],df[i,'FirstName'],'Certificate',sep="_")
  out_filename = paste0(out_filename, ifelse(word_out, '.docx','.pdf'))

  #save customized Rmd to a temporary file
  write_file(current_cert, "tmp.Rmd")

  #create the certificates using R markdown.
  #it will detect the ending of the output file and use the right format
  rmarkdown::render("tmp.Rmd", output_file = out_filename)

  #temporary Rmd file can be deleted.
  file.remove("tmp.Rmd")

}



