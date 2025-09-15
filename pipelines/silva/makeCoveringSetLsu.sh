#kmerfilterset.sh in=lsu_deduped100pct.fa.gz k=31 rcomp=f out=lsu_covering_31mers.fa maxkpp=1
#reformat.sh in=lsu_deduped100pct.fa.gz out=lsu_deduped100pct_padded.fa.gz padleft=31 padright=31 ow
#shred.sh in=lsu_deduped100pct_padded.fa.gz length=150 minlength=62 overlap=30 out=shreds.fa.gz ow
#kmerfilterset.sh in=shreds.fa.gz initial=lsu_covering_31mers.fa k=31 rcomp=f out=lsu_shred_covering_31mers.fa maxkpp=1
nohup time kmerfilterset.sh in=shreds.fa.gz initial=lsu_shred_covering_31mers_temp.fa k=31 rcomp=f out=lsu_shred_covering_31mers.fa maxkpp=1 maxpasses=40000 fastawrap=99999 -Xmx2g ow
