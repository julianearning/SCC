plot(df[1,], type="l", ylim=c(-0.1,0.1))
df <- as.matrix(df)

plot(df[1,], type="l", ylim=c(-0.1,0.1))


for(i in seq(1,nrow(df), 100)) {
	plot(df[i,], type="b", ylim=c(-0.1,0.1))
	Sys.sleep(0.7)
}
