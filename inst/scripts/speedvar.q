provide.data(aircraft)
hw <- nnbr(Speed, 30)
hw <- hw/exp(mean(log(hw)))
par(mfrow=c(1,2))
sm.density(Speed, yht=0.0022, positive=T)
sm.density(Speed, yht=0.0022, xlim=c(-700,4000), h.weights=hw)
par(mfrow=c(1,1))

