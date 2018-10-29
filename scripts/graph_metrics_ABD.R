library(igraph)
library(ggplot2)
library(lme4)
library(arm)

subj  <- as.factor(c('0050952', '0050953', '0050954', '0050955', '0050956', '0050957', '0050958', '0050959', '0050960', '0050961', '0050962', '0050964', '0050965', '0050966', '0050967', '0050968', '0050969', '0050970', '0050971', '0050972', '0050973', '0050974', '0050975', '0050976', '0050977', '0050978', '0050979', '0050980', '0050981', '0050982', '0050983', '0050984', '0050985', '0050986', '0050987', '0050988', '0050989', '0050990', '0050991', '0050992', '0050993', '0050994', '0050995', '0050996', '0050997', '0050998', '0050999', '0051000', '0051001', '0051002', '0051003', '0051006', '0051007', '0051008', '0051009', '0051010', '0051011', '0051012', '0051013', '0051014', '0051015', '0051016', '0051017', '0051018', '0051019', '0051020', '0051021', '0051023', '0051024', '0051025', '0051026', '0051027', '0051028', '0051029', '0051030', '0051032', '0051033', '0051034', '0051035', '0051036', '0051038', '0051039', '0051040', '0051041', '0051042', '0051044', '0051045', '0051046', '0051047', '0051048', '0051049', '0051050', '0051051', '0051052', '0051053', '0051054', '0051055', '0051056', '0051057', '0051058', '0051059', '0051060', '0051061', '0051062', '0051063', '0051064', '0051065', '0051066', '0051067', '0051068', '0051069', '0051070', '0051072', '0051073', '0051074', '0051075', '0051076', '0051077', '0051078', '0051079', '0051080', '0051081', '0051082', '0051083', '0051084', '0051085', '0051086', '0051087', '0051088', '0051089', '0051090', '0051091', '0051093', '0051094', '0051095', '0051096', '0051097', '0051098', '0051099', '0051100', '0051101', '0051102', '0051103', '0051104', '0051105', '0051106', '0051107', '0051108', '0051109', '0051110', '0051111', '0051112', '0051113', '0051114', '0051115', '0051116', '0051117', '0051118', '0051119', '0051120', '0051121', '0051122', '0051123', '0051124', '0051125', '0051126', '0051127', '0051128', '0051129', '0051130', '0051131', '0051146', '0051147', '0051148', '0051149', '0051150', '0051151', '0051152', '0051153', '0051154', '0051155', '0051156', '0051159'))


### FUNCTIONS ###

# participation coefficient
# (adjacency, community)
part.coef <- function(w,Ci) {
	n <- nrow(w) # number of rows (and columns as matrix should be symmetrical 
	Ko <- apply(w,1,sum) # grab node-wise degree
	pos.val <- w # copy matrix to mask
	pos.val[pos.val!=0] <- 1 # binary matrix for community mask
	Gc <- pos.val%*%diag(c(Ci)) # create community mask
	Kc2 <- vector('numeric',n) # container for community-wise degree
	for (i in 1:max(Ci)) {
		Kc2 = Kc2+rowSums(w*(Gc==i))^2 # grab community degree
	}
	P <- round(array(1,n)-Kc2/(Ko^2),4) # calculate and round participation coefficient
	return(P)
}

# read individual adjacency matrices
# (subject vector, directory, number of regions, condition)
read.ind.mat <- function(sub,dir,n) { 
	setwd(dir)
	fin <- array(0, dim=c(n,n,length(sub)))
	ii <- 0 # index
	for (s in levels(sub)) {
		cat(s,'\n')
		ii <- ii+1
		temp <- as.matrix(read.table(paste0(dir,'/',s,'YEO.adj.txt'))) # read matrix
		diag(temp) <- 0 # set diagonal to 0
		temp[temp<0] <- 0 # set negatives to 0
		fin[,,ii] <- temp
	}
	return(fin)
}

# Group averaged matrices
grp.mat <- function(sub,dir,n) {
	setwd(dir)
	fin <- array(0, dim=c(n,n))
	ii <- 0 # index
	for (s in levels(sub)) {
		cat(s,'\n')
		ii <- ii+1
		temp <- as.matrix(read.table(paste0(dir,'/',s,'YEO.adj.txt'))) # read matrix
		diag(temp) <- 0 # set diagonal to 0
		fin <- fin+temp
	}
	fin <- fin/length(sub)
	return(fin)
}

### SETUP ###

dir <- '~/Dropbox/DSCN/Experiments/ABIDE/network_analysis/imfar_2017/yeo_adj_mat'
# Ci <- as.matrix(read.table('~/iCloud/CHT/scripts/power_communities.txt'))
Ci <- read.csv('~/Dropbox/DSCN/Experiments/ABIDE/network_analysis/imfar_2017/yeo_vol_com.csv')
pheno <- read.csv('~/Dropbox/DSCN/Experiments/ABIDE/network_analysis/imfar_2017/NYU_pheno.csv',colClasses='character')
pheno$Mean <- as.numeric(pheno$Mean)
pheno$SRS_RAW_TOTAL <- as.numeric(pheno$SRS_RAW_TOTAL)
pheno$FIQ <- as.numeric(pheno$FIQ)
pheno$ADOS_TOTAL <- as.numeric(pheno$ADOS_TOTAL)
pheno$ADOS_SOCIAL <- as.numeric(pheno$ADOS_SOCIAL)
pheno$AGE_AT_SCAN <- as.numeric(pheno$AGE_AT_SCAN)
ind<- read.ind.mat(subj,dir,128)
graph <- list()

for (ss in 1:length(subj)) {
	graph[[subj[ss]]] <- graph_from_adjacency_matrix(ind[,,ss], mode='undirected', weighted=TRUE)
}

Q <- vector('numeric',length(subj))
for (ss in 1:length(subj)) {
	Q[ss] <- modularity(graph[[subj[ss]]],Ci$com)
}
Q <- data.frame(Q=Q,SUB_ID=subj)

fin <- merge(pheno,Q,by='SUB_ID')

P <- array(0, dim=c(128,length(subj)))
for (ss in 1:length(subj)) {
	P[,ss] <- part.coef(ind[,,ss],Ci$com)
}
P <- aperm(P)
P <- data.frame(P,SUB_ID=subj)


roi <- 'X108'
tP <- data.frame(SUB_ID=P$SUB_ID,roi=P[,roi])
new <- merge(fin,tP,by='SUB_ID')
good <- subset(new, new$SEX==1 & new$Mean<0.1 & new$AGE_AT_SCAN<19.74)
good$SRS_RAW_TOTAL[good$SRS_RAW_TOTAL==-9999] <- NA
mod <- lm(roi ~ DX_GROUP*AGE_AT_SCAN + Mean + FIQ + EYE_STATUS_AT_SCAN, data=good)
summary(mod)

mod <- lm(roi ~ Mean + FIQ + EYE_STATUS_AT_SCAN, data=good)
new <- data.frame(res=residuals(mod),grp=good$DX_GROUP,age=good$AGE_AT_SCAN)
ggplot(new, aes(age,res,color=grp)) + geom_point() + geom_smooth(method='lm')


#### IMFAR 2017 figures ####
> a <- delete.edges(a, which(E(a)$weight < threshold)-1)

ggplot(new, aes(age,res,color=grp)) + 
	geom_point() + 
	geom_smooth(method='lm',fill='grey80') +
	scale_color_manual(values=c('forestgreen','steelblue')) +
	scale_x_continuous(breaks=c(6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
	scale_y_continuous(breaks=c(-.1,0)) +
	theme_bw() +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text.x=element_text(size=14),
		axis.text.y=element_text(size=14))


asd <- subset(good, good$DX_GROUP==1)
td <- subset(good, good$DX_GROUP==2)
mod2 <- lm(roi ~ AGE_AT_SCAN + Mean + FIQ + EYE_STATUS_AT_SCAN + ADOS_SOCIAL, data=asd)
summary(mod2)

ggplot(good, aes(AGE_AT_SCAN,Q,color=DX_GROUP)) + geom_point() + geom_smooth(method='lm')
good$DX_GROUP <- factor(good$DX_GROUP, levels=c(2,1))
mod <- lm(Q ~ DX_GROUP + AGE_AT_SCAN + Mean + FIQ + EYE_STATUS_AT_SCAN, data=good)
summary(mod)
s <- sim(mod,10000)
ef <- data.frame(s@coef)

#### IMFAR 2017 figures #### 

ggplot(good, aes(DX_GROUP,Q,fill=DX_GROUP)) + 
	geom_jitter(width=.4) +
	geom_boxplot(alpha=.75) +
	scale_fill_manual(values=c('steelblue','forestgreen')) +
	theme_bw() +
	labs(x='',y='') +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.text.y=element_text(size=16))



ggplot(ef, aes(DX_GROUP1)) + 
	geom_vline(xintercept=0,size=1) +
	geom_density(alpha=0.6, fill = 'darkorchid') +
	labs(x='',y='') +
	theme_bw() +
	theme(legend.position='none',
		axis.text.x=element_text(size=16),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank())






mod <- lm(Q ~ AGE_AT_SCAN + Mean + FIQ + EYE_STATUS_AT_SCAN, data=good)
new <- data.frame(res=residuals(mod),grp=good$DX_GROUP)
ggplot(new, aes(grp,res)) + geom_boxplot() + geom_jitter(width=0.6)

# X117 - rTPJ
# X108 - lTPJ


#### group maps ####

m <- mean(good$AGE_AT_SCAN)
young <- subset(good,good$AGE_AT_SCAN<=m)
old <- subset(good,good$AGE_AT_SCAN>m)
roi <- c(108) # 117

a.y <- subset(young,young$DX_GROUP==1)
gASD.y <- grp.mat(as.factor(a.y$SUB_ID),dir,128)
gASD.y[-roi,-roi] <- 0

a.o <- subset(old,old$DX_GROUP==1)
gASD.o <- grp.mat(as.factor(a.o$SUB_ID),dir,128)
gASD.o[-roi,-roi] <- 0

t.y <- subset(young,young$DX_GROUP==2)
gTD.y <- grp.mat(as.factor(t.y$SUB_ID),dir,128)
gTD.y[-roi,-roi] <- 0

t.o <- subset(old,old$DX_GROUP==2)
gTD.o <- grp.mat(as.factor(t.o$SUB_ID),dir,128)
gTD.o[-roi,-roi] <- 0

guide <- read.table('~/iCloud/ABD/cir_guide.txt',stringsAsFactors=FALSE)
names(guide) <- c('name','on','off')

# prepare links for circos
prep.circos.links <- function(adj,guide) {
	temp <- adj*lower.tri(adj)
	temp <- round(temp,3)
	N <- ncol(adj)
	n <- length(temp[temp>0])
	s.name <- vector('character',n)
	s.on <- vector('numeric',n)
	s.off <- vector('numeric',n)
	e.name <- vector('numeric',n)
	e.on <- vector('numeric',n)
	e.off <- vector('numeric',n)
	type <- vector('character',n)
	ind <- 1
	for (ii in 1:N) {
		cat('...',ii)
		for (jj in 1:N) {
			if (temp[ii,jj]>0) {
				s.name[ind] <- guide$name[ii]
				s.on[ind] <- guide$on[ii]
				s.off[ind] <- guide$off[ii]
				e.name[ind] <- guide$name[jj]
				e.on[ind] <- guide$on[jj]
				e.off[ind] <- guide$off[jj]
				type[ind] <- paste0('type=1,score=',temp[ii,jj])
				ind <- ind+1
			}
		}
	}
	fin <- data.frame(s.name,s.on,s.off,e.name,e.on,e.off,type)
	return(fin)
}

# links <- prep.circos.links(gASD.y,guide)
# write.table(links,'~/iCloud/ABD/circos/asdY/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# links <- prep.circos.links(gASD.o,guide)
# write.table(links,'~/iCloud/ABD/circos/asdO/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# links <- prep.circos.links(gTD.y,guide)
# write.table(links,'~/iCloud/ABD/circos/tdY/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# links <- prep.circos.links(gTD.o,guide)
# write.table(links,'~/iCloud/ABD/circos/tdO/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)


# td.asd.y <- gTD.y-gASD.y
# links <- prep.circos.links(td.asd.y,guide)
# write.table(links,'~/iCloud/ABD/circos/td-asdY/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# td.asd.o <- gTD.o-gASD.o
# links <- prep.circos.links(td.asd.o,guide)
# write.table(links,'~/iCloud/ABD/circos/td-asdO/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# asd.td.y <- gASD.y-gTD.y
# links <- prep.circos.links(asd.td.y,guide)
# write.table(links,'~/iCloud/ABD/circos/asd-tdY/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)

# asd.td.o <- gASD.o-gTD.o
# links <- prep.circos.links(asd.td.o,guide)
# write.table(links,'~/iCloud/ABD/circos/asd-tdO/data/links.txt',row.names=FALSE,col.names=FALSE,quote=FALSE)


#### Test participation in left TPJ in younger half on sample
young <- subset(good, good$AGE_AT_SCAN < 12.8)
old <- subset(good, good$AGE_AT_SCAN > 12.8)
t.test(roi ~ DX_GROUP,data=young)


############ Threshold
t.1 <- ind
t.1[t.1<0.1] <-0
t.2 <- ind
t.2[t.2<0.2] <- 0
t.3 <- ind
t.3[t.3<0.3] <- 0
t.4 <- ind
t.4[t.4<0.4] <- 0
t.5 <- ind
t.5[t.5<0.5] <- 0

P.1 <- array(0, dim=c(128,length(subj)))
P.2 <- array(0, dim=c(128,length(subj)))
P.3 <- array(0, dim=c(128,length(subj)))
P.4 <- array(0, dim=c(128,length(subj)))
P.5 <- array(0, dim=c(128,length(subj)))

for (ss in 1:length(subj)) {
	P.1[,ss] <- part.coef(t.1[,,ss],Ci$com)
}
P.1 <- aperm(P.1)
P.1 <- data.frame(P.1,SUB_ID=subj)

for (ss in 1:length(subj)) {
	P.2[,ss] <- part.coef(t.2[,,ss],Ci$com)
}
P.2 <- aperm(P.2)
P.2 <- data.frame(P.2,SUB_ID=subj)

for (ss in 1:length(subj)) {
	P.3[,ss] <- part.coef(t.3[,,ss],Ci$com)
}
P.3 <- aperm(P.3)
P.3 <- data.frame(P.3,SUB_ID=subj)

for (ss in 1:length(subj)) {
	P.4[,ss] <- part.coef(t.4[,,ss],Ci$com)
}
P.4 <- aperm(P.4)
P.4 <- data.frame(P.4,SUB_ID=subj)

for (ss in 1:length(subj)) {
	P.5[,ss] <- part.coef(t.5[,,ss],Ci$com)
}
P.5 <- aperm(P.5)
P.5 <- data.frame(P.5,SUB_ID=subj)

t <- data.frame(p=P[,108],p1=P.1[,108],p2=P.2[,108],p3=P.3[,108],p4=P.4[,108],p5=P.5[,108],SUB_ID=new$SUB_ID,DX_GROUP=new$DX_GROUP,AGE_AT_SCAN=new$AGE_AT_SCAN,Mean=new$Mean,EYE_STATUS_AT_SCAN=new$EYE_STATUS_AT_SCAN,FIQ=new$FIQ,SEX=new$SEX)
good <- subset(t, t$SEX==1 & t$Mean<0.1 & t$AGE_AT_SCAN<19.74)
mod <- lm(p3 ~ DX_GROUP*AGE_AT_SCAN + Mean + FIQ + EYE_STATUS_AT_SCAN, data=good)
summary(mod)



