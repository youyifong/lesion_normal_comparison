a="CD3"
counts = mysapply(slides, function (s) 
  c(
    cnt_total=nrow(subset(alldat[[s]][[a]])),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<500 & epi &    region=="intact")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<500 & dermis & region=="intact")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<500 & epi &    region=="ulcerative")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<500 & dermis & region=="ulcerative")),
    cnt_total_500=nrow(subset(alldat[[s]][[a]], dist2DEJ<500)),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<1000 & epi &    region=="intact")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<1000 & dermis & region=="intact")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<1000 & epi &    region=="ulcerative")),
    nrow(subset(alldat[[s]][[a]], dist2DEJ<1000 & dermis & region=="ulcerative")),
    cnt_total_1000=nrow(subset(alldat[[s]][[a]], dist2DEJ<1000))
  )
)
counts=as.data.frame(counts)
counts$initial = get.initial(slides)
counts$shedding_rate = severity$Shedding.rate[match(counts$initial, severity$initial)]
counts

# scatterplot of shedding_rate vs cnt_total
png("output/figures/correlation_shedding_rate_CD3.png", width=10, height=10, res=300)


par(mfrow=c(1,3))
reduce_margin()
plot(counts$shedding_rate, counts$cnt_total, 
     xlab="Shedding rate", ylab="Total cells", main=a, ylim=c(0,15000))
# put R2 in legend
r2 = summary(lm(counts$cnt_total ~ counts$shedding_rate))$r.squared
legend("topright", legend=paste("R2=", round(r2, 3)), bty="n")
# abline(lm(counts$cnt_total ~ counts$shedding_rate), col="red")

plot(counts$shedding_rate, counts$cnt_total_1000, 
     xlab="Shedding rate", ylab="Total cells (dist2DEJ<1000)", main=a, ylim=c(0,15000))
r2 = summary(lm(counts$cnt_total_1000 ~ counts$shedding_rate))$r.squared
legend("topright", legend=paste("R2=", round(r2, 3)), bty="n")
# abline(lm(counts$cnt_total_1000 ~ counts$shedding_rate), col="red")

plot(counts$shedding_rate, counts$cnt_total_500, 
     xlab="Shedding rate", ylab="Total cells (dist2DEJ<500)", main=a, ylim=c(0,15000))
r2 = summary(lm(counts$cnt_total_500 ~ counts$shedding_rate))$r.squared
legend("topright", legend=paste("R2=", round(r2, 3)), bty="n")
# abline(lm(counts$cnt_total_500 ~ counts$shedding_rate), col="red")
