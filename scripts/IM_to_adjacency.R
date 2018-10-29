library(Hmisc)
subj <- as.factor(c("RED_LL_167", "RED_LL_187", "RED_RL_106", "RED_RL_111", "RED_RL_118", "RED_RL_119", "RED_RL_124", "RED_RL_126", "RED_RL_128", "RED_RL_130", "RED_RL_146", "RED_RL_150", "RED_RL_151", "RED_RL_152", "RED_RL_153", "RED_RL_154", "RED_RL_155", "RED_RL_157", "RED_RL_158", "RED_RL_159", "RED_RL_160", "RED_RL_162", "RED_RL_163", "RED_RL_165"))
con <- as.factor(c("peer","comp"))
dir <- "~/iCloud/CHT/betas"
out_dir <- "~/iCloud/CHT/sig_adj_mat"

for (s in levels(subj)){
	cat(s,'\n')
	for (c in levels(con)){
		cat('  ',c,'\n')
		raw <- read.table(paste0(dir,"/",s,".",c,".stats.1D"), stringsAsFactors=FALSE)[-1,-1]
		raw <- matrix(unlist(lapply(raw,as.numeric)), ncol = length(raw), byrow=FALSE)
		pear.mat <- rcorr(raw)
		fin.mat <- (pear.mat$P<0.01)*(pear.mat$r)
		fin.mat[fin.mat<0] <- 0
		fin.mat <- round(fin.mat, 4)
		write.table(fin.mat, paste0(out_dir,"/",s,".",c,".thresh.adj.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
	}
}
