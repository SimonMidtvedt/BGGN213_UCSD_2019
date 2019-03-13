Lecture 7. Functions and Packages
================
Simon Midtvedt
1/30/2019

### Funcions revisited

Load (i.e. **source**) our rescale() funcion from last day.

``` r
source("http://tinyurl.com/rescale-R")
```

Test this function

``` r
rescale(1:5)
```

    ## [1] 0.00 0.25 0.50 0.75 1.00

``` r
#rescale(c(1:5, "string"))
```

Which would give an error. We want to make this function more acceptable.

``` r
rescale2(c(1:5, "string"))
```

Which gives a describing error message.

Now we want to make a function that finds NAs on the same spot in two vectors. We will start with forgetting about the function and focus on the task with an example.

``` r
x <- c(1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
```

Use Google to find a built-in function that may help

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

``` r
which(is.na(x) & is.na(y))
```

    ## [1] 3

Now we take the working snippet and make a first function

``` r
both_na <- function(x, y) {
  # Check for NA elements in both input vectors
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x, y)
```

    ## [1] 1

``` r
x <-  c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
y3 <- c(1, NA, NA, NA, NA)
```

``` r
# What will this return?
both_na(x, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
#both_na2(x, y)
```

Gives the following message: Error: Input x and y should be vectors of the same length.

``` r
both_na3(x, y)
```
