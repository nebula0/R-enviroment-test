a <- 40

test_func <- function(){
    #assign("a", 30, envir = .GlobalEnv)
    a <<- 20
}

test_func()

a
