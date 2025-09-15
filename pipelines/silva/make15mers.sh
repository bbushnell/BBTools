#23S:
#Make initial set
kmercountexact.sh in=lsu_deduped100pct.fa.gz mincount=5000 out=lsu_mincount5000_k15.fa k=15 rcomp=f ow

#Test coverage
bbduk.sh in=lsu_deduped100pct.fa.gz ref=lsu_mincount5000_k15.fa k=15 mm=f rcomp=f

#Make other set
kmercountexact.sh in=ssu_deduped100pct.fa.gz mincount=10 out=ssu_mincount10_k15.fa k=15 rcomp=f ow

#Subtract sets
bbduk.sh ref=ssu_mincount10_k15.fa in=lsu_mincount5000_k15.fa k=15 mm=f outm=shared_10_k15.fa outu=lsuonly_10_k15.fa ow

#Test coverage again
bbduk.sh in=lsu_deduped100pct.fa.gz ref=lsuonly_10_k15.fa k=15 mm=f rcomp=f

#Compress
kcompress.sh in=lsuonly_10_k15.fa k=15 out=lsu_15mers.fa.gz ow zl=11 rcomp=f

#Test coverage one last time
bbduk.sh in=lsu_deduped100pct.fa.gz ref=lsu_15mers.fa.gz k=15 mm=f rcomp=f



#16S:
#Make initial set
kmercountexact.sh in=ssu_deduped100pct.fa.gz mincount=5000 out=ssu_mincount5000_k15.fa k=15 rcomp=f ow

#Test coverage
bbduk.sh in=ssu_deduped100pct.fa.gz ref=ssu_mincount5000_k15.fa k=15 mm=f rcomp=f

#Make other set
kmercountexact.sh in=lsu_deduped100pct.fa.gz mincount=8 out=lsu_mincount8_k15.fa k=15 rcomp=f ow

#Subtract sets
bbduk.sh ref=lsu_mincount8_k15.fa in=ssu_mincount5000_k15.fa k=15 mm=f outm=shared_8_k15.fa outu=ssuonly_8_k15.fa ow

#Test coverage again
bbduk.sh in=ssu_deduped100pct.fa.gz ref=ssuonly_8_k15.fa k=15 mm=f rcomp=f

#Compress
kcompress.sh in=ssuonly_8_k15.fa k=15 out=ssu_15mers.fa.gz ow zl=11 rcomp=f

#Test coverage one last time
bbduk.sh in=ssu_deduped100pct.fa.gz ref=ssu_15mers.fa.gz k=15 mm=f rcomp=f



#5S: TODO


#Reverse test
bbduk.sh in=lsu_deduped100pct.fa.gz ref=ssu_15mers.fa.gz k=15 mm=f rcomp=f

#Reverse test
bbduk.sh in=ssu_deduped100pct.fa.gz ref=lsu_15mers.fa.gz k=15 mm=f rcomp=f
