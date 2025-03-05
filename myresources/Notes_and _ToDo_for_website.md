# Todo

* Update this document to reflect move to Quarto.
* Fix logos on project page

# Notes 

## General
* Some information that applies to this site and handelgroup is not repeated here, see the handelgroup notes file.

* It seems that qmd sites can't handle the weaving of code with markdown using the read_chunk setup. They ignore eval=FALSE and run stuff anyway. Might need to run those as Rmd files until fixed. - SEEMS FIXED NOW?

* Using these Quarto extensions: fontawsome, include-code-files


## Needed R packages

We use renv to manage the packages needed. It usually finds automatically what is needed. As a backup, this is a hopefully complete list of all R packages needed to recompile all the posts of the website:

install.packages(c('cowplot','geosphere','scholar','wordcloud','bibliometrix','tidytext','visdat','kableExtra','janitor','ggplot2','dplyr','stringr','tidytuesdayR','readr','emoji',"rethinking","cmdstanr","brms","cmdstanr"))

remotes::install_github('rmcelreath/rethinking')
remotes::install_github('andreashandel/flowdiagramr')

`cmdstanr` and `cmdstan` are iffy and often misbehave. 


## Website Logo
* Save your icon as a square 512x512 pixel image named icon.png and place the image in your root assets/images/ folder, creating the assets and images folders if they don't already exist. (https://sourcethemes.com/academic/docs/customization/)
* To show logo on site, place icon/logo in assets/images and/or assets/media


## Redirects
Follow the instructions in my Hugo -> Quarto conversion blog post using the post-render script setup.

## To build website
Run: quarto render in console

## To deploy site
Run this quarto build/publish command on console:

quarto publish netlify

This should automatically try to connect to Netlify and publish.

No automated building from Github at this point, need to rebuild with quarto, then run quarto command to push updates to Netlify.


## To implement commenting
I'm using utterances, implemented by installing utterances into github folder and adding comments.html into the layouts/partials folder
https://mscipio.github.io/post/utterances-comment-engine/
https://masalmon.eu/2019/10/02/disqus/
https://www.davidfong.info/post/hugoacademiccommentswithutterances/



## Change footer



## Ideas for new posts/content



### Science/Research

* Reading/managing/publishing papers: turn my presentation and my text on the 'teaching-research-resources' list into a blog post.

* Blog post discussing stats vs mechanistic models and non-parametric vs parametric (how they are somewhat similar, one trades off power for assumptions). Almost like a bias variance trade-off.

* How to approach your first (and every subsequent) research project – use text I wrote for ‘onboarding’. Also combine with research project slides/talk.

*	Post Writing a scientific article with R/Quarto/etc.

*	Turn some presentations, e.g. ‘good research project’ into blog posts.

### Teaching/Tech

Review of my technology stack for online courses
  * Github
  * Slack
  * gradingapp

* Experiences teaching a course using Github organization versus not.
  
* gh_class package to manage class with github and my experience using github for teaching: workflow, exercises

* conceptual thought of usefulness of certain teaching aids like metaphors/examples/jokes that might confuse-enlighten or disengage-engage

* quizgrader description


### Others

* See "started blog posts": staying mentally in good shape in grad school (and in general). comparing mental exercise/meditation to physical exercise, where do analogies work and where do they not.

* epidemiological biases and "learn from experts" literature (e.g. tim ferriss books): Use/consume/learn from folks in Tim Ferriss blog and similar, while knowing the huge impact of survivor bias and non-random sampling. (Example: Jim Collins wisdom on Tim Ferriss show and his likely poor analysis (The Model Thinker, p 162)

* Fully implement, then write about book review system: Excel - rmarkdown/blogdown

* Blog on ‘not reconciling differences’, e.g., Use/insist on the best (statistical) methodology possible, while knowing/appreciating that it often doesn’t matter.

* Impact/usefulness of different activities (research/teaching/outreach/mentoring, etc.)