for (name in c("small","moderate")) {
 s <- if(name=="small") rep(25L,4) else rep(100L,20); m<-sum(s); packed<-sum(s*(s+1)/2)
 print(data.frame(workload=name,n=if(name=="small")200L else 2000L,m=m,blocks=length(s),
  min_block=min(s),median_block=median(s),max_block=max(s),
  filter=c("hard_truncate","ridge_fixed","ridge_lw"),tau=.01,eta=.1,
  stored_bytes=packed*4+2*m*4+m*8,
  largest_transient_block_bytes=max(s)^2*8*4+max(s)*(if(name=="small")200L else 2000L)*4,
  sweep_value_visits=sum(s^2),rebuild_value_visits=sum(s^2)),row.names=FALSE)
}
cat("Analytical regression signal only; no peak-RSS or performance claim.\n")
