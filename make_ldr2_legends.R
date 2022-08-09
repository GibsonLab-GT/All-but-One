library(readr)

# Make color legends for LD R2 values per peak

pdf(file = "legends.pdf")
# make: peak 1
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#00008B", "#0000BE", "#1616F9", "#4443ED", "#9F9CD5"))
mtext("LD R2 Peak 1", at=0.2, cex=2)
# make: peak 2
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#006400", "#4CA814", "#99ED29", "#B7ED62", "#C5D5A6"))
mtext("LD R2 Peak 2", at=0.2, cex=2)
# make: peak 3
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#CD950C", "#E3C046", "#F9EB80", "#EEE7A2", "#E3DDAF"))
mtext("LD R2 Peak 3", at=0.2, cex=2)
# make: peak 4
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#8B1A1A", "#BE3A2E", "#D84A38", "#F96E55", "#E39B8F"))
mtext("LD R2 Peak 4", at=0.2, cex=2)
# make: peak 5
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#551A8B", "#9040AC", "#CB66CD", "#D58DD1", "#CFB5CB"))
mtext("LD R2 Peak 5", at=0.2, cex=2)
# make: peak 6
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#8B4513", "#BE6F2D", "#F29A48", "#EEB177", "#D8C1AD"))
mtext("LD R2 Peak 6", at=0.2, cex=2)
# make: peak 7
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#8B0A50", "#BE5A86", "#F2ABBD", "#EEC3CA", "#D8C7C9"))
mtext("LD R2 Peak 7", at=0.2, cex=2)
# make: peak 8
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#8B4500", "#BE6F00", "#F29A00", "#EEB143","#D8C19C"))
mtext("LD R2 Peak 8", at=0.2, cex=2)
# make: peak 9
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#458B74", "#5EBE9E", "#78F2C9", "#99EDD0", "#BBD5CB"))
mtext("LD R2 Peak 9", at=0.2, cex=2)
# make: peak 10
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("1.0", "0.8", "0.6", "0.4", "0.2"), pch=16, pt.cex=3, cex=1.5, bty='n',
       col = c("#473C8B", "#6253BE", "#7D6AF2", "#9C8DED", "#BCB5D5"))
mtext("LD R2 Misc", at=0.2, cex=2)
dev.off()
