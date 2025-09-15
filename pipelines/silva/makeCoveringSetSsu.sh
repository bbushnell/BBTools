kmerfilterset.sh in=ssu_deduped100pct.fa.gz k=31 rcomp=f out=ssu_covering_31mers.fa maxkpp=1 -Xmx8g
reformat.sh in=ssu_deduped100pct.fa.gz out=ssu_deduped100pct_padded.fa.gz padleft=31 padright=31 ow
shred.sh in=ssu_deduped100pct_padded.fa.gz length=150 minlength=62 overlap=30 out=shredsSsu.fa.gz ow
time kmerfilterset.sh in=shredsSsu.fa.gz initial=ssu_covering_31mers.fa k=31 rcomp=f out=ssu_shred_covering_31mers.fa maxkpp=1 maxpasses=200000 fastawrap=99999 -Xmx8g ow
