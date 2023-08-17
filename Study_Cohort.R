load("/home/liuwent/04-Full_Model/pd.RData")

PE_cases <- pd[pd$Sample_Group=='Disease',]
dim(PE_cases)
# 24 14

controls <- pd[pd$Sample_Group=='Controls',]
dim(controls)
# 38 14

## Age:
mean(PE_cases$Age)
sd(PE_cases$Age)

mean(controls$Age)
sd(controls$Age)

## Parity:
mean(PE_cases$Parity)
sd(PE_cases$Parity)

mean(controls$Parity)
sd(controls$Parity)

## BMI:
mean(PE_cases$BMI)
sd(PE_cases$BMI)

mean(controls$BMI)
sd(controls$BMI)

## GA:
mean(PE_cases$GA)
sd(PE_cases$GA)

mean(controls$GA)
sd(controls$GA)

## Smoker:
table(PE_cases$Smoker)
table(controls$Smoker)

## Eth2:
table(PE_cases$Eth2)
table(controls$Eth2)

# t-test:
t.test(PE_cases$Age, controls$Age)
t.test(PE_cases$Parity, controls$Parity)
t.test(PE_cases$BMI, controls$BMI)
t.test(PE_cases$GA, controls$GA)

# chisq.test:
chisq.test(table(PE_cases$Smoker), table(controls$Smoker))








