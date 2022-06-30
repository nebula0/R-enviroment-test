library(openxlsx)
df <- read.xlsx("t_test.xlsx")


for(i in 1:10) {
    print(i)
    print(names(df[i]))
    print(shapiro.test(df[,i]))
}
# p抓0.01可以的:地毯草、兩耳草

var.test(df[,"A地毯草"], df[,"B地毯草"])
# 地毯草變異數同
var.test(df[,"A兩耳草"], df[,"B兩耳草"])
# 
df[,"A兩耳草"]
A_two_ear <- c(31,  76,  83,   3,  39,   8,  14, 105,  26,  11,  20,  10)
df[,"B兩耳草"]
B_two_ear <- c(40,  5, 14, 13,  0,  0, 16,  9,  3, 16,  5, 16)
var.test(A_two_ear, B_two_ear)
# 兩耳草變異數不同 p<0.01
A_carpet <- df[,"A地毯草"]
B_carpet <- df[,"B地毯草"]
df[,"A地毯草"]
A_carpet
var.test(A_carpet, B_carpet)
# 地毯草變異數相同 p>0.01

t.test(A_two_ear, B_two_ear, var.equal = FALSE)
# p>0.01無顯著差異

t.test(A_carpet, B_carpet, var.equal = TRUE)
# p<0.01 有顯著差異

