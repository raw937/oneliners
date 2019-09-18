library(reshape2)
library(ggplot2)

phyla <- read.delim("taxa_phyla_krona.txt", check.names=F,stringsAsFactors=F)
phage <- read.delim("taxa_phage-phyla_krona.txt", check.names=F,stringsAsFactors=F)

mm <- melt(phyla)
mm1 <- melt(phage)

#ggplot2 barplot
c <- ggplot(mm, aes(x=reorder(phyla,-value), value, fill=variable)) 
c <- c + geom_bar(stat="identity", position=position_dodge())
c <- c + theme_bw(base_size=15, base_family="helvetica") #removes background/sets up text size
c <- c + theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
c <- c + scale_fill_brewer(palette = "Set2")
c <- c + labs(x="", y="Relative abundance (%)", title="")
c <- c + guides(fill=guide_legend(title=NULL))
c <- c + theme(legend.justification=c(1,0), legend.position=c(1,.75))
c

#ggplot2 barplot
c1 <- ggplot(mm1, aes(x=reorder(phyla,-value), value, fill=variable)) 
c1 <- c1 + geom_bar(stat="identity", position=position_dodge())
c1 <- c1 + theme_bw(base_size=15, base_family="helvetica") #removes background/sets up text size
c1 <- c1 + theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
c1 <- c1 + scale_fill_brewer(palette = "Set2")
c1 <- c1 + labs(x="", y="Relative abundance (%)", title="")
c1 <- c1 + guides(fill=guide_legend(title=NULL))
c1 <- c1 + theme(legend.justification=c(1,0), legend.position=c(1,.75))
c1


#ggplot2 on its side
#ggplot2 barplot
c <- ggplot(mm, aes(x=reorder(phyla,+value), value, fill=variable)) 
c <- c + geom_bar(stat="identity", position=position_dodge())
c <- c + theme_bw(base_size=15, base_family="helvetica") #removes background/sets up text size
c <- c + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
c <- c + scale_fill_brewer(palette = "Set2")
c <- c + labs(x="", y="Relative abundance (%)", title="")
c <- c + guides(fill=guide_legend(title=NULL))
c <- c + theme(legend.justification=c(1,0), legend.position=c(1,.75))
c <- c + coord_flip()
c

#ggplot2 barplot
c1 <- ggplot(mm1, aes(x=reorder(phyla,+value), value, fill=variable)) 
c1 <- c1 + geom_bar(stat="identity", position=position_dodge())
c1 <- c1 + theme_bw(base_size=15, base_family="helvetica") #removes background/sets up text size
c1 <- c1 + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
c1 <- c1 + scale_fill_brewer(palette = "Set2")
c1 <- c1 + labs(x="", y="Relative abundance (%)", title="")
c1 <- c1 + guides(fill=guide_legend(title=NULL))
c1 <- c1 + theme(legend.justification=c(1,0), legend.position=c(1,.1))
c1 <- c1 + coord_flip()
c1


#pie chart

df <- data.frame(
  group = c("Lysogenic", "Lytic"),
  value = c(64, 36)
)
head(df)

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie