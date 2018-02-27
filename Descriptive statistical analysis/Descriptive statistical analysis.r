library(fBasics)
da <- read.table("CNY_EUR.txt",header = T)
dim(da)
ced <- da[,3]
basicStats(ced)

##对收盘价计算对数收益率并查看其趋势
r<- numeric(length(ced)-1)
for(i in 2:length(ced)){
  r[i-1] <- (log(ced[i])-log(ced[i-1]))
}
basicStats(r)
#转化为时间序列
rts <- ts(r)
plot(rts,xlab="",ylab="对数收益率")

#绘制对数收益率直方图
hist(r,breaks = 100,main = "收益率分布图",xlab = "收益率",ylab = "频数")

#绘制对数收益率的QQ图
qqnorm(r,main = "QQ图",xlab = "理论分位数",ylab = "样本分位数", plot.it = TRUE)
qqline(r)

#Jarque-Bera检验
jarqueberaTest(r)

#自相关性和偏自相关性
acf(r,main = "")
pacf(r,main = "")

#Ljung-Box检验
Box.test(r,lag = 10, type = "Ljung-Box")

#收益率序列平方自相关性和偏自相关性
acf(r*r,main = "")
pacf(r*r,main = "")

#对收益率平方序列进行Ljung-Box检验
Box.test(r*r,lag = 10, type = "Ljung-Box")


