plot(0, type="n", ylab="", xlab="")
segments(.7, .7,  1.3, .7, col="red", lwd=4)
segments(.7, .5,  1.3, .5, col="blue", lwd=4)
segments(1, .0,  1.3, .0, col="red", lwd=4)
segments(.7, .0,  1., .0, col="blue", lwd=4)

segments(1, -.5,  1.3, -.5, col="red", lwd=4)
segments(.7, -.5,  1., -.5, col="blue", lwd=4)
segments(.7, -.5,  .8, -.5, col="red", lwd=4)

text(.65, .7, labels = "maternal", cex=.5)
text(.65, .5, labels = "paternal", cex=.5)
text(.65, .0, labels = "single\nrecombinant", cex=.5)
text(.65, -.5, labels = "double\nrecombinant", cex=.5)