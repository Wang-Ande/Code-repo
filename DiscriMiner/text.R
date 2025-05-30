library(DiscriMiner)
# Input Normalization Data
# Data format:samples in rows and variables in columns.

# data input
data <- read.csv("./01_Data/01.MetQuant/meta_intensity_combined.csv",row.names = 1)
data <- data[,c(14:40)]
data <- log2(data)
data <- t(data)

# group input
group <- read.xlsx("./01_Data/01.MetQuant/sam_infor_combined.xlsx")
group <- group[c(1:27),]

# PLS-DA
target_data <- data[grep("MOLM13",rownames(data)),]
target_group <- group[grep("MOLM13",group$id),3]


my_pls1 <- plsDA(variables = target_data, group = target_group, 
                 autosel = FALSE, comps = 2,validation = NULL, learn = NULL, 
                 test = NULL,cv = "LOO", k = NULL, retain.models = FALSE)
plot(my_pls1)

my_pls1$classification
my_pls1$error_rate
my_pls1$components
my_pls1$Q2
my_pls1$R2
my_pls1$VIP
str(my_pls1)


library(ggplot2)
q2_values <- data.frame(Component = rownames(my_pls1$Q2), Q2 = my_pls1$Q2[, "Q2.global"])
r2_values <- data.frame(Component = rownames(my_pls1$R2), R2 = my_pls1$R2[, "R2Ycum"])

ggplot(q2_values, aes(x = Component, y = Q2)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = round(Q2, 2)), vjust = -0.3) +
  labs(title = "Q2 Values", x = "Component", y = "Q2") +
  theme_minimal()

ggplot(r2_values, aes(x = Component, y = R2)) +
  geom_bar(stat = "identity", fill = "green") +
  geom_text(aes(label = round(R2, 2)), vjust = -0.3) +
  labs(title = "R2 Values", x = "Component", y = "R2") +
  theme_minimal()

