library(data.table)

header <- c("UID:604591897", "email:rsumbaly@ucla.edu", "Undergrad or Grad:Grad")

# read tab separated SNP data
SNP_data <- read.table("SNP_Status.txt", header = TRUE)

# separate SNP data into cases and controls
cases <- subset(SNP_data, Status == 'Case')
controls <- subset(SNP_data, Status == 'Control')

# get number of cases, controls and SNPs
n_total <- dim(SNP_data)[1]
n_cases <- dim(cases)[1]
n_controls <- dim(controls)[1]
n_SNPs <- dim(SNP_data)[2] - 1

# alpha, bonferroni correction calculation 
alpha <- 0.05
bonferroni_corr <- alpha / n_SNPs

# prepare output variables
A <- c("<A>")
B <- c("<B>")
C <- c("<C>")

SNP_name <- colnames(SNP_data)
chi_squared <- numeric()
p_values <- numeric()
SA <- numeric()

# calculate association statistic for each SNP 
# apply bonferroni correction to obtain significant SNPs
for (SNP in 1:n_SNPs) {
  # obtain p+, p- and p for each SNP
  p_plus <- sum(cases[SNP])/(2 * n_cases)
  p_minus <- sum(controls[SNP])/(2 * n_controls)
  p_avg <- (p_plus + p_minus) / 2
  
  association_statistic <- (p_plus - p_minus)/(sqrt(2/n_total) * sqrt(p_avg * (1 - p_avg)))
  SA <- c(SA, association_statistic)
  
  if (!is.na(association_statistic)) {
    chi_squared <- c(chi_squared, association_statistic**2)  
  }
  
  if (p_avg == 0) {
    p_value <- 1.00
  }
  else {
    p_value <- 2 * pnorm(abs(association_statistic), lower.tail = FALSE) # two-tailed   
  }
  p_values <- c(p_values, p_value)
  A <- c(A, paste(SNP_name[SNP], p_value, sep = ":"))

  if (p_value < bonferroni_corr) {
    B <- c(B, SNP_name[SNP])
  }
}

# calculating inflation statistics
df <- 1
median_df <- df * ( 1 - (2 / (9 * df))**3 )
lambda_gc <- median(chi_squared)/ median_df

A <- c(A, "</A>")
B <- c(B, "</B>")
C <- c(C, paste("Lambda_gc", lambda_gc, sep = ":"), "</C>")

# Q-Q Plots
# Association Statistics
qqplot(qnorm(ppoints(length(SA))), SA, main = expression("Q-Q plot for Standard Normal"), 
             xlab = "Expected Values - SA", ylab = "Observed Values - SA")
qqline(SA)

# Chi-Squared Statistic
qqplot(qchisq(ppoints(length(chi_squared)), df = 1), y = chi_squared, main = expression("Q-Q plot for" ~~ {chi^2}[nu == 1]) , 
              xlab = "Expected Values - Chi-Squared Statistic", ylab = "Observed Values - Chi-Squared Statistic")
qqline(chi_squared)

# output to file
output_file <- "PA1_Output.txt"

# remove output file if exists
if(file.exists(output_file))
  file.remove(output_file)

writeLines(c(header,A,B,C), file(output_file))
