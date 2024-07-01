# Based on: https://blog.djnavarro.net/posts/2022-04-20_porting-to-quarto/#netlify-redirects
# lines to insert to a netlify _redirect file
redirects <- paste0("/talk /presentations \n/talks /presentations \n/presentations/ /presentations \n/post /posts")
# write the _redirect file
writeLines(redirects, here::here("_site", "_redirects"))