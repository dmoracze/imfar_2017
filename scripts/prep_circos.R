str <- read.csv('~/iCloud/ABD/cir_prep_str.csv',stringsAsFactors=FALSE)
#str$module <- factor(str$module,levels=c('Sub','VisCent','VisPeri','SomMotA','SomMotB','LimbicA','LimbicB','ContA','ContB','ContC','VentAttnA','VentAttnB','DorsAttnA','DorsAttnB','DefaultA','DefaultB','DefaultC','DefaultD'))
str$module <- factor(str$module,levels=c('Subcortical','Visual','Somatomotor','Limbic','Control','VentAttn','DorsAttn','Default'))
N <- length(str$module)
hemi <- as.factor(c('lh','rh'))

# structure.label.txt
chr <- vector('character',n)
init <- vector('numeric',n)
end <- vector('numeric',n)
name <- vector('character',n)

ii <- 0
for (mm in levels(str$module)) {
	l.loc <- 0
	r.loc <- 0
	for (hh in levels(hemi)) {
		temp <- subset(str,str$module==mm & str$hemi==hh)
		for (rr in 1:length(temp$roi)) {
			ii <- ii+1
			if (temp$hemi[rr]=='lh') {
				chr[ii] <- paste0(temp$module[rr],'-l')
				init[ii] <- l.loc
				end[ii] <- l.loc+99
				name[ii] <- temp$roi[rr]
				l.loc <- l.loc+100

			} else {
				chr[ii] <- paste0(temp$module[rr],'-r')
				init[ii] <- r.loc
				end[ii] <- r.loc+99
				name[ii] <- temp$roi[rr]
				r.loc <- r.loc+100
			}
		}
	}
}
str.label <- data.frame(chr,init,end,name)
write.table(str.label,'~/iCloud/ABD/asd/data/structure.label.txt',col.names=FALSE,row.names=FALSE,quote=FALSE)

# segments
dat <- str.label
col.1 <- vector('character',n)
col.2 <- vector('character',n)
col.3 <- vector('character',n)
col.4 <- vector('character',n)
col.5 <- vector('numeric',n)
col.6 <- vector('numeric',n)
col.7 <- vector('character',n)

ii <- 0
for (mm in levels(dat$chr)) {
	ii <- ii+1
	temp <- subset(dat, dat$chr==mm)
	col.1[ii] <- 'chr'
	col.2[ii] <- '-'
	col.3[ii] <- as.character(temp$chr[1])
	col.4[ii] <- strsplit(as.character(temp$chr[1]),'-')[[1]][1]
	col.5[ii] <- 0
	col.6[ii] <- max(temp$end)
	col.7[ii] <- 'black'
	for (rr in 1:length(temp$name)) {
		ii <- ii+1
		col.1[ii] <- 'band'
		col.2[ii] <- as.character(temp$chr[rr])
		col.3[ii] <- as.character(temp$name[rr])
		col.4[ii] <- as.character(temp$name[rr])
		col.5[ii] <- temp$init[rr]
		col.6[ii] <- temp$end[rr]
		col.7[ii] <- tolower(temp$name[rr])
	}
}
seg <- data.frame(col.1,col.2,col.3,col.4,col.5,col.6,col.7)
write.table(seg,'~/iCloud/ABD/asd/data/segments.txt',col.names=FALSE,row.names=FALSE,quote=FALSE)

# color.brain.conf


write.table(col,'~/iCloud/ABD/asd/etc/color.brain.conf',col.names=FALSE,row.names=FALSE,quote=FALSE)

