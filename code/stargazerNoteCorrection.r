x <- capture.output(stargazer:::.stargazer.wrap)
idx <- c(grep("for \\(i in 1:length\\(\\.format\\.cutoffs\\)\\)", x)[2],
         grep("if \\(!is\\.null\\(notes\\)\\)", x),
         grep("if \\(!is\\.null\\(notes\\.align\\)\\)", x)[2])
eval(parse(text = paste0("stargazer.wrap <- ", paste(x[c(1:(idx[1] - 1), 
                                                         (idx[2]):(idx[3] - 1), 
                                                         idx[1]:(idx[2] - 1), 
                                                         idx[3]:(length(x) - 2))], collapse = "\n"))))

## Create a new stargazer.() that uses our modified stargazer.wrap() function
x <- capture.output(stargazer)
x <- gsub(".stargazer.wrap", "stargazer.wrap", x)
eval(parse(text = paste0("stargazer. <- ", paste(x[1:(length(x)-2)], collapse = "\n"))))
