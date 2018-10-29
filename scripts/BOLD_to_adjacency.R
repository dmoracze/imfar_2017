library(Hmisc)

subj  <- as.factor(c('0050952', '0050953', '0050954', '0050955', '0050956', '0050957', '0050958', '0050959', '0050960', '0050961', '0050962', '0050964', '0050965', '0050966', '0050967', '0050968', '0050969', '0050970', '0050971', '0050972', '0050973', '0050974', '0050975', '0050976', '0050977', '0050978', '0050979', '0050980', '0050981', '0050982', '0050983', '0050984', '0050985', '0050986', '0050987', '0050988', '0050989', '0050990', '0050991', '0050992', '0050993', '0050994', '0050995', '0050996', '0050997', '0050998', '0050999', '0051000', '0051001', '0051002', '0051003', '0051006', '0051007', '0051008', '0051009', '0051010', '0051011', '0051012', '0051013', '0051014', '0051015', '0051016', '0051017', '0051018', '0051019', '0051020', '0051021', '0051023', '0051024', '0051025', '0051026', '0051027', '0051028', '0051029', '0051030', '0051032', '0051033', '0051034', '0051035', '0051036', '0051038', '0051039', '0051040', '0051041', '0051042', '0051044', '0051045', '0051046', '0051047', '0051048', '0051049', '0051050', '0051051', '0051052', '0051053', '0051054', '0051055', '0051056', '0051057', '0051058', '0051059', '0051060', '0051061', '0051062', '0051063', '0051064', '0051065', '0051066', '0051067', '0051068', '0051069', '0051070', '0051071', '0051072', '0051073', '0051074', '0051075', '0051076', '0051077', '0051078', '0051079', '0051080', '0051081', '0051082', '0051083', '0051084', '0051085', '0051086', '0051087', '0051088', '0051089', '0051090', '0051091', '0051093', '0051094', '0051095', '0051096', '0051097', '0051098', '0051099', '0051100', '0051101', '0051102', '0051103', '0051104', '0051105', '0051106', '0051107', '0051108', '0051109', '0051110', '0051111', '0051112', '0051113', '0051114', '0051115', '0051116', '0051117', '0051118', '0051119', '0051120', '0051121', '0051122', '0051123', '0051124', '0051125', '0051126', '0051127', '0051128', '0051129', '0051130', '0051131', '0051146', '0051147', '0051148', '0051149', '0051150', '0051151', '0051152', '0051153', '0051154', '0051155', '0051156', '0051159'))



#subj <- as.factor(c('0051071'))
index <- read.csv('~/Dropbox/DSCN/Experiments/ABIDE/network_analysis/reordered_abide_rois_SURF.csv')

dat_dir <- c('~/iCloud/ABD/NYU')
out_dir <- c('~/iCloud/ABD/yeo_adj_mat')

for (s in levels(subj)){
	cat(paste(s,'\n'))
	# ventral striatum
	dim <- read.table(paste(dat_dir,'/vol/',s,'.VS.1D', sep=''),stringsAsFactors=F)[-1,-1]
	# right amygdala
	amy.rh <- as.numeric(read.table(paste(dat_dir,'/vol/',s,'.rh.amy.1D', sep=''),stringsAsFactors=F)[-1,-1])
	# left amygdala
	amy.lh <- as.numeric(read.table(paste(dat_dir,'/vol/',s,'.lh.amy.1D', sep=''),stringsAsFactors=F)[-1,-1])
	# # yeo 17 right
	yeo.rh <- read.table(paste(dat_dir,'/yeo/',s,'_rh.yeo17split.1D',sep=''),stringsAsFactors=F)[-1,-1]
	# # yeo 17 left
	yeo.lh <- read.table(paste(dat_dir,'/yeo/',s,'_lh.yeo17split.1D', sep=''),stringsAsFactors=F)[-1,-1]
	# # assemble into timeseries matrix
	raw <- cbind(dim,amy.rh,amy.lh,yeo.rh,yeo.lh)
	raw <- matrix(unlist(lapply(raw,as.numeric)),ncol=length(raw),byrow=F)
	nraw <- raw[,index$new.reind]
	# create adjacency matrix
	cormat <- rcorr(raw)$r
	cormat <- round(cormat,4)
	diag(cormat) <- 0
	cormat[cormat<0] <- 0
	write.table(cormat, paste(out_dir,'/',s,'YEO.adj.txt',sep=''),row.names=F,col.names=F,quote = F)
}

# for (s in levels(subj)){
# 	cat(paste(s,'\n'))
# 	# # yeo 17 right
# 	yeo.rh <- read.table(paste(dat_dir,'/yeo/',s,'_rh.yeo17split.1D',sep=''),stringsAsFactors=F)[-1,-1]
# 	# # yeo 17 left
# 	yeo.lh <- read.table(paste(dat_dir,'/yeo/',s,'_lh.yeo17split.1D', sep=''),stringsAsFactors=F)[-1,-1]
# 	# # assemble into timeseries matrix
# 	raw <- cbind(yeo.rh,yeo.lh)
# 	raw <- matrix(unlist(lapply(raw,as.numeric)),ncol=length(raw),byrow=F)
# 	nraw <- raw[,index$new.reind]
# 	# create adjacency matrix
# 	cormat <- rcorr(raw)$r
# 	cormat <- round(cormat,4)
# 	write.table(cormat, paste(out_dir,'/',s,'.surf.adj.txt',sep=''),row.names=F,col.names=F,quote = F)
# }











