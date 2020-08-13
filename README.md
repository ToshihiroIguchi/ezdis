## ezdis
Easy distribution fitting

### Description
GUI univariate distribution fitting and comparison.

### Launch web application
First, access [CRAN](https://cran.r-project.org/), download R and install.
You can install from R console.
First install related packages.

    install.packages(c("shiny", "tibble", "dplyr", "DT", "readr", "readxl", "R.utils", "fitdistrplus", "ggplot2", "ismev", "evd", "actuar", "FAdist", "extraDistr", "EnvStats", "RcppFaddeeva", "mixtools", "goftest"))
    
If packages are installed, it can be started from R console with the following command.
    
    shiny::runGitHub("ezdis", "ToshihiroIguchi")

### Host
Host the Shiny application from GitHub in a private network.
Enter the following command in R console.

    #Port specification
    port <- 6500

    #Acquire private address information
    ipconfig.dat <- system("ipconfig", intern = TRUE)
    ipv4.dat <- ipconfig.dat[grep("IPv4", ipconfig.dat)][1]
    ip <- gsub(".*? ([[:digit:]])", "\\1", ipv4.dat)

    #Host the Shiny application from GitHub
    shiny::runGitHub("ezdis", "ToshihiroIguchi", launch.browser = FALSE, port = port, host = ip)

If you are in the private network, you can also launch the Shiny application by entering the URL following `Listing on` to the browser.

### shinyapps.io
Hosted on shinyapps.io.
Use time is limited.
https://toshihiroiguchi.shinyapps.io/ezdis/

### How to use
`ezdis` can be operated on the browser.
You can load csv, xls, xlsx and xlsm format files by clicking the `Browse...` button at the top left of the screen.

The first line is the item name, and only the columns with all numerical values are read.
However, if there is a sheet of xls, xlsx, xlsm format and all have numerical values, all the numerical values of that sheet are analyzed.

This time, analyze the `sample.csv` file.
After uploading to ezdis, you can select `Use data`.
By selecting the item name, the survey will start to find out which distribution fits.
`ezdis` investigates whether it fits the distribution selected in `Setting` and displays it in `Distribution` in ascending order of AIC.


### License 

```
MIT License

Copyright (c) 2020 Toshihiro Iguchi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

### Auther
Toshihiro Iguchi