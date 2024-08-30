printTitle <- function(x, fontcase = c("upper", "lower", "proper")) {
     fontcase <- match.arg(fontcase)
     x <- ifelse(fontcase == "proper", toproper(x), ifelse(fontcase == "upper", toupper(x), tolower(x)))
     cat("\n\n", x, "\n", sep = "")
     cat(rep("-", nchar(x)), sep = "")
     cat("\n")
}


toproper <- function(x) {
     s <- strsplit(x, " ")[[1]]
     paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
           sep = "", collapse = " ")
}
