library(RegCalReliab)

############# EX #################

main = read.table("/Users/liubowen/Desktop/EPI515_BOWEN_REVISED/lecture3/data/main_part2.dat", header = T)
reliab = read.table("/Users/liubowen/Desktop/EPI515_BOWEN_REVISED/lecture3/data/ers_part2.dat", header = T)


#----------------Logistic_EX--------------------

ex_log = RC_ExReliab(formula = comb ~ sbp(sbp2, sbp3),
                     main_data = main,
                     rep_data = reliab,
                     link = "logistic")

ex_log$uncorrected
ex_log$corrected



#----------------Linear_EX--------------------

ex_lin = RC_ExReliab(formula = dbp ~ sbp(sbp2, sbp3) + chol(chol2,chol3),
                     main_data = main,
                     rep_data = reliab,
                     link = "linear")

ex_lin$uncorrected
ex_lin$corrected


#----------------Poisson_EX--------------------

ex_poisson = RC_ExReliab(formula = comb ~ sbp(sbp2, sbp3) + chol(chol2, chol3) + gluc(gluc2, gluc3) + bmi(bmi2, bmi3) +
                       smoke + ages2 + ages3 + ages4,
                     main_data = main,
                     rep_data = reliab,
                     link = "log")

ex_poisson$uncorrected
ex_poisson$corrected


############# IN #################
library(RegCalReliab)
main_in = read.table("/Users/liubowen/Desktop/EPI515_BOWEN_REVISED/lecture3/data/internal_data.dat", header = T)
main_in[is.na(main_in)] = 0


#----------------Logistic_IN--------------------


in_log = RC_InReliab(formula = comb ~ mysbp(sbp,sbp2, sbp3) + chol,
                     main_data = main_in,
                     link = "logistic")

in_log$uncorrected
in_log$corrected


#----------------Linear_IN--------------------

in_lin = RC_InReliab(formula = dbp ~ mysbp(sbp,sbp2, sbp3) + mychol(chol,chol2, chol3),
                     main_data = main_in,
                     link = "linear")

in_lin$uncorrected
in_lin$corrected


#----------------Poisson_IN--------------------

in_poisson = RC_InReliab(formula = comb ~ mysbp(sbp,sbp2, sbp3) + chol,
                     main_data = main_in,
                     link = "log")

in_poisson$uncorrected
in_poisson$corrected




