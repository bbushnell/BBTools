#Link a copy of deduplicated Sliva (LSU and SSU trunc combined) with Us changed to Ts.
ln -s /global/projectb/sandbox/gaag/bbtools/silva/latest/both_deduped_sorted.fa.gz .

#Aggresively remove similar sequences
dedupe.sh in=both_deduped_sorted.fa.gz out=dd2_e2_s5.fa.gz s=5 e=2 ordered zl=9 fastawrap=4000

#Less-aggressively remove duplicate sequences
clumpify.sh in=both_deduped_sorted.fa.gz out=clumped_s2.fa.gz zl=9 s=2 dedupe fastawrap=4000 ow passes=2

#Entropy-mask the sequences
bbduk.sh -Xmx1g in=dd2_e2_s5.fa.gz out=dd2_e2_s5_masked.fa.gz zl=9 entropy=0.6 entropyk=4 entropywindow=24 maskentropy ordered ow qtrim=rl trimq=1 fastawrap=4000
bbduk.sh -Xmx1g in=clumped_s2.fa.gz out=clumped_s2_masked.fa.gz zl=9 entropy=0.6 entropyk=4 entropywindow=24 maskentropy ordered ow qtrim=rl trimq=1 fastawrap=4000

#Generate synthetic reads
randomreads.sh -Xmx31g adderrors=f ref=clumped_s2_masked.fa.gz reads=300m out=synth_s2.fa.gz len=100 zl=6 fastawrap=4000 illuminanames
randomreads.sh -Xmx8g adderrors=f ref=dd2_e2_s5_masked.fa.gz reads=300m out=synth_e2_s5.fa.gz len=100 zl=6 fastawrap=4000 illuminanames

#Remove duplicate reads
clumpify.sh in=synth_s2.fa.gz out=synth_s2_clumped_s1.fa.gz reorder zl=9 fastawrap=4000 groups=1 dedupe s=1 rcomp
clumpify.sh in=synth_e2_s5.fa.gz out=synth_e2_s5_clumped_s1.fa.gz reorder zl=9 fastawrap=4000 groups=1 dedupe s=1 rcomp

#Create baseline kmer sets at different depths for different sensitivites (only one depth is needed)
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=1000 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers1000A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=500 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers500A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=200 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers200A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=100 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers100A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=50 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers50A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=40 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers40A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=30 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers30A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=20 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers20A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=10 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers10A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=8 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers8A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=5 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers5A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=3 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers3A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=2 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers2A.fa.gz
kcompress.sh -Xmx31g ow zl=9 pigz=16 min=1 in=dd2_e2_s5_masked.fa.gz out=stdout.fa | clumpify.sh -Xmx16g in=stdin.fa k=16 reorder groups=1 fastawrap=4000 ow zl=9 pigz=16 out=riboKmers1A.fa.gz

#Find the missed synthetic reads
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers1000A.fa.gz out=read_misses1000A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers500A.fa.gz out=read_misses500A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers200A.fa.gz out=read_misses200A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers100A.fa.gz out=read_misses100A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers50A.fa.gz out=read_misses50A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers40A.fa.gz out=read_misses40A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers30A.fa.gz out=read_misses30A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers20A.fa.gz out=read_misses20A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers10A.fa.gz out=read_misses10A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers8A.fa.gz out=read_misses8A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers5A.fa.gz out=read_misses5A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers3A.fa.gz out=read_misses3A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers2A.fa.gz out=read_misses2A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
bbduk.sh -Xmx8g in=synth_e2_s5_clumped_s1.fa.gz ref=riboKmers1A.fa.gz out=read_misses1A.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow

#Iterate over each depth to add missed kmers (again, only one depth is needed)
kcompress.sh -Xmx31g ow min=2000 in=read_misses1000A.fa.gz out=riboKmers1000B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses1000A.fa.gz ref=riboKmers1000A.fa.gz,riboKmers1000B.fa.gz out=read_misses1000B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=1000 in=read_misses1000B.fa.gz out=riboKmers1000C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers1000A.fa.gz,riboKmers1000B.fa.gz,riboKmers1000C.fa.gz out=read_misses1000C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=500 in=read_misses1000C.fa.gz out=riboKmers1000D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses1000C.fa.gz ref=riboKmers1000A.fa.gz,riboKmers1000B.fa.gz,riboKmers1000C.fa.gz,riboKmers1000D.fa.gz out=read_misses1000D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=250 in=read_misses1000D.fa.gz out=riboKmers1000E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses1000D.fa.gz ref=riboKmers1000A.fa.gz,riboKmers1000B.fa.gz,riboKmers1000C.fa.gz,riboKmers1000D.fa.gz,riboKmers1000E.fa.gz out=read_misses1000E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers1000A.fa.gz,riboKmers1000B.fa.gz,riboKmers1000C.fa.gz,riboKmers1000D.fa.gz,riboKmers1000E.fa.gz out=riboKmers1000merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers1000merged.fa.gz out=riboKmers1000clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers1000clumped.fa.gz out=riboKmers1000fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=1000 in=read_misses500A.fa.gz out=riboKmers500B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses500A.fa.gz ref=riboKmers500A.fa.gz,riboKmers500B.fa.gz out=read_misses500B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=500 in=read_misses500B.fa.gz out=riboKmers500C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers500A.fa.gz,riboKmers500B.fa.gz,riboKmers500C.fa.gz out=read_misses500C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=250 in=read_misses500C.fa.gz out=riboKmers500D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses500C.fa.gz ref=riboKmers500A.fa.gz,riboKmers500B.fa.gz,riboKmers500C.fa.gz,riboKmers500D.fa.gz out=read_misses500D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=125 in=read_misses500D.fa.gz out=riboKmers500E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses500D.fa.gz ref=riboKmers500A.fa.gz,riboKmers500B.fa.gz,riboKmers500C.fa.gz,riboKmers500D.fa.gz,riboKmers500E.fa.gz out=read_misses500E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers500A.fa.gz,riboKmers500B.fa.gz,riboKmers500C.fa.gz,riboKmers500D.fa.gz,riboKmers500E.fa.gz out=riboKmers500merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers500merged.fa.gz out=riboKmers500clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers500clumped.fa.gz out=riboKmers500fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=500 in=read_misses200A.fa.gz out=riboKmers200B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses200A.fa.gz ref=riboKmers200A.fa.gz,riboKmers200B.fa.gz out=read_misses200B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=250 in=read_misses200B.fa.gz out=riboKmers200C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers200A.fa.gz,riboKmers200B.fa.gz,riboKmers200C.fa.gz out=read_misses200C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=125 in=read_misses200C.fa.gz out=riboKmers200D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses200C.fa.gz ref=riboKmers200A.fa.gz,riboKmers200B.fa.gz,riboKmers200C.fa.gz,riboKmers200D.fa.gz out=read_misses200D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=60 in=read_misses200D.fa.gz out=riboKmers200E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses200D.fa.gz ref=riboKmers200A.fa.gz,riboKmers200B.fa.gz,riboKmers200C.fa.gz,riboKmers200D.fa.gz,riboKmers200E.fa.gz out=read_misses200E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers200A.fa.gz,riboKmers200B.fa.gz,riboKmers200C.fa.gz,riboKmers200D.fa.gz,riboKmers200E.fa.gz out=riboKmers200merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers200merged.fa.gz out=riboKmers200clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers200clumped.fa.gz out=riboKmers200fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=200 in=read_misses100A.fa.gz out=riboKmers100B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses100A.fa.gz ref=riboKmers100A.fa.gz,riboKmers100B.fa.gz out=read_misses100B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=100 in=read_misses100B.fa.gz out=riboKmers100C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers100A.fa.gz,riboKmers100B.fa.gz,riboKmers100C.fa.gz out=read_misses100C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=50 in=read_misses100C.fa.gz out=riboKmers100D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses100C.fa.gz ref=riboKmers100A.fa.gz,riboKmers100B.fa.gz,riboKmers100C.fa.gz,riboKmers100D.fa.gz out=read_misses100D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=25 in=read_misses100D.fa.gz out=riboKmers100E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses100D.fa.gz ref=riboKmers100A.fa.gz,riboKmers100B.fa.gz,riboKmers100C.fa.gz,riboKmers100D.fa.gz,riboKmers100E.fa.gz out=read_misses100E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers100A.fa.gz,riboKmers100B.fa.gz,riboKmers100C.fa.gz,riboKmers100D.fa.gz,riboKmers100E.fa.gz out=riboKmers100merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers100merged.fa.gz out=riboKmers100clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers100clumped.fa.gz out=riboKmers100fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=100 in=read_misses50A.fa.gz out=riboKmers50B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses50A.fa.gz ref=riboKmers50A.fa.gz,riboKmers50B.fa.gz out=read_misses50B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=50 in=read_misses50B.fa.gz out=riboKmers50C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers50A.fa.gz,riboKmers50B.fa.gz,riboKmers50C.fa.gz out=read_misses50C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=25 in=read_misses50C.fa.gz out=riboKmers50D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses50C.fa.gz ref=riboKmers50A.fa.gz,riboKmers50B.fa.gz,riboKmers50C.fa.gz,riboKmers50D.fa.gz out=read_misses50D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=12 in=read_misses50D.fa.gz out=riboKmers50E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses50D.fa.gz ref=riboKmers50A.fa.gz,riboKmers50B.fa.gz,riboKmers50C.fa.gz,riboKmers50D.fa.gz,riboKmers50E.fa.gz out=read_misses50E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers50A.fa.gz,riboKmers50B.fa.gz,riboKmers50C.fa.gz,riboKmers50D.fa.gz,riboKmers50E.fa.gz out=riboKmers50merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers50merged.fa.gz out=riboKmers50clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers50clumped.fa.gz out=riboKmers50fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=80 in=read_misses40A.fa.gz out=riboKmers40B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses40A.fa.gz ref=riboKmers40A.fa.gz,riboKmers40B.fa.gz out=read_misses40B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=40 in=read_misses40B.fa.gz out=riboKmers40C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers40A.fa.gz,riboKmers40B.fa.gz,riboKmers40C.fa.gz out=read_misses40C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=20 in=read_misses40C.fa.gz out=riboKmers40D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses40C.fa.gz ref=riboKmers40A.fa.gz,riboKmers40B.fa.gz,riboKmers40C.fa.gz,riboKmers40D.fa.gz out=read_misses40D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=10 in=read_misses40D.fa.gz out=riboKmers40E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses40D.fa.gz ref=riboKmers40A.fa.gz,riboKmers40B.fa.gz,riboKmers40C.fa.gz,riboKmers40D.fa.gz,riboKmers40E.fa.gz out=read_misses40E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers40A.fa.gz,riboKmers40B.fa.gz,riboKmers40C.fa.gz,riboKmers40D.fa.gz,riboKmers40E.fa.gz out=riboKmers40merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers40merged.fa.gz out=riboKmers40clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers40clumped.fa.gz out=riboKmers40fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=60 in=read_misses30A.fa.gz out=riboKmers30B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses30A.fa.gz ref=riboKmers30A.fa.gz,riboKmers30B.fa.gz out=read_misses30B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=30 in=read_misses30B.fa.gz out=riboKmers30C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers30A.fa.gz,riboKmers30B.fa.gz,riboKmers30C.fa.gz out=read_misses30C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=15 in=read_misses30C.fa.gz out=riboKmers30D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses30C.fa.gz ref=riboKmers30A.fa.gz,riboKmers30B.fa.gz,riboKmers30C.fa.gz,riboKmers30D.fa.gz out=read_misses30D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=8 in=read_misses30D.fa.gz out=riboKmers30E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses30D.fa.gz ref=riboKmers30A.fa.gz,riboKmers30B.fa.gz,riboKmers30C.fa.gz,riboKmers30D.fa.gz,riboKmers30E.fa.gz out=read_misses30E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers30A.fa.gz,riboKmers30B.fa.gz,riboKmers30C.fa.gz,riboKmers30D.fa.gz,riboKmers30E.fa.gz out=riboKmers30merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers30merged.fa.gz out=riboKmers30clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers30clumped.fa.gz out=riboKmers30fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=40 in=read_misses20A.fa.gz out=riboKmers20B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses20A.fa.gz ref=riboKmers20A.fa.gz,riboKmers20B.fa.gz out=read_misses20B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=20 in=read_misses20B.fa.gz out=riboKmers20C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers20A.fa.gz,riboKmers20B.fa.gz,riboKmers20C.fa.gz out=read_misses20C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=12 in=read_misses20C.fa.gz out=riboKmers20D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses20C.fa.gz ref=riboKmers20A.fa.gz,riboKmers20B.fa.gz,riboKmers20C.fa.gz,riboKmers20D.fa.gz out=read_misses20D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=7 in=read_misses20D.fa.gz out=riboKmers20E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses20D.fa.gz ref=riboKmers20A.fa.gz,riboKmers20B.fa.gz,riboKmers20C.fa.gz,riboKmers20D.fa.gz,riboKmers20E.fa.gz out=read_misses20E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers20A.fa.gz,riboKmers20B.fa.gz,riboKmers20C.fa.gz,riboKmers20D.fa.gz,riboKmers20E.fa.gz out=riboKmers20merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers20merged.fa.gz out=riboKmers20clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers20clumped.fa.gz out=riboKmers20fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=20 in=read_misses10A.fa.gz out=riboKmers10B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses10A.fa.gz ref=riboKmers10A.fa.gz,riboKmers10B.fa.gz out=read_misses10B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=12 in=read_misses10B.fa.gz out=riboKmers10C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers10A.fa.gz,riboKmers10B.fa.gz,riboKmers10C.fa.gz out=read_misses10C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=8 in=read_misses10C.fa.gz out=riboKmers10D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses10C.fa.gz ref=riboKmers10A.fa.gz,riboKmers10B.fa.gz,riboKmers10C.fa.gz,riboKmers10D.fa.gz out=read_misses10D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=5 in=read_misses10D.fa.gz out=riboKmers10E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses10D.fa.gz ref=riboKmers10A.fa.gz,riboKmers10B.fa.gz,riboKmers10C.fa.gz,riboKmers10D.fa.gz,riboKmers10E.fa.gz out=read_misses10E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers10A.fa.gz,riboKmers10B.fa.gz,riboKmers10C.fa.gz,riboKmers10D.fa.gz,riboKmers10E.fa.gz out=riboKmers10merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers10merged.fa.gz out=riboKmers10clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers10clumped.fa.gz out=riboKmers10fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=16 in=read_misses8A.fa.gz out=riboKmers8B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses8A.fa.gz ref=riboKmers8A.fa.gz,riboKmers8B.fa.gz out=read_misses8B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=10 in=read_misses8B.fa.gz out=riboKmers8C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers8A.fa.gz,riboKmers8B.fa.gz,riboKmers8C.fa.gz out=read_misses8C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=7 in=read_misses8C.fa.gz out=riboKmers8D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses8C.fa.gz ref=riboKmers8A.fa.gz,riboKmers8B.fa.gz,riboKmers8C.fa.gz,riboKmers8D.fa.gz out=read_misses8D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=4 in=read_misses8D.fa.gz out=riboKmers8E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses8D.fa.gz ref=riboKmers8A.fa.gz,riboKmers8B.fa.gz,riboKmers8C.fa.gz,riboKmers8D.fa.gz,riboKmers8E.fa.gz out=read_misses8E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers8A.fa.gz,riboKmers8B.fa.gz,riboKmers8C.fa.gz,riboKmers8D.fa.gz,riboKmers8E.fa.gz out=riboKmers8merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers8merged.fa.gz out=riboKmers8clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers8clumped.fa.gz out=riboKmers8fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=12 in=read_misses5A.fa.gz out=riboKmers5B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses5A.fa.gz ref=riboKmers5A.fa.gz,riboKmers5B.fa.gz out=read_misses5B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=8 in=read_misses5B.fa.gz out=riboKmers5C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers5A.fa.gz,riboKmers5B.fa.gz,riboKmers5C.fa.gz out=read_misses5C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=5 in=read_misses5C.fa.gz out=riboKmers5D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses5C.fa.gz ref=riboKmers5A.fa.gz,riboKmers5B.fa.gz,riboKmers5C.fa.gz,riboKmers5D.fa.gz out=read_misses5D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=3 in=read_misses5D.fa.gz out=riboKmers5E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses5D.fa.gz ref=riboKmers5A.fa.gz,riboKmers5B.fa.gz,riboKmers5C.fa.gz,riboKmers5D.fa.gz,riboKmers5E.fa.gz out=read_misses5E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers5A.fa.gz,riboKmers5B.fa.gz,riboKmers5C.fa.gz,riboKmers5D.fa.gz,riboKmers5E.fa.gz out=riboKmers5merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers5merged.fa.gz out=riboKmers5clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers5clumped.fa.gz out=riboKmers5fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=10 in=read_misses3A.fa.gz out=riboKmers3B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses3A.fa.gz ref=riboKmers3A.fa.gz,riboKmers3B.fa.gz out=read_misses3B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=7 in=read_misses3B.fa.gz out=riboKmers3C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers3A.fa.gz,riboKmers3B.fa.gz,riboKmers3C.fa.gz out=read_misses3C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=4 in=read_misses3C.fa.gz out=riboKmers3D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses3C.fa.gz ref=riboKmers3A.fa.gz,riboKmers3B.fa.gz,riboKmers3C.fa.gz,riboKmers3D.fa.gz out=read_misses3D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=3 in=read_misses3D.fa.gz out=riboKmers3E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses3D.fa.gz ref=riboKmers3A.fa.gz,riboKmers3B.fa.gz,riboKmers3C.fa.gz,riboKmers3D.fa.gz,riboKmers3E.fa.gz out=read_misses3E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers3A.fa.gz,riboKmers3B.fa.gz,riboKmers3C.fa.gz,riboKmers3D.fa.gz,riboKmers3E.fa.gz out=riboKmers3merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers3merged.fa.gz out=riboKmers3clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers3clumped.fa.gz out=riboKmers3fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

kcompress.sh -Xmx31g ow min=8 in=read_misses2A.fa.gz out=riboKmers2B.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses2A.fa.gz ref=riboKmers2A.fa.gz,riboKmers2B.fa.gz out=read_misses2B.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=6 in=read_misses2B.fa.gz out=riboKmers2C.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=synth_s2_clumped_s1.fa.gz ref=riboKmers2A.fa.gz,riboKmers2B.fa.gz,riboKmers2C.fa.gz out=read_misses2C.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=4 in=read_misses2C.fa.gz out=riboKmers2D.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses2C.fa.gz ref=riboKmers2A.fa.gz,riboKmers2B.fa.gz,riboKmers2C.fa.gz,riboKmers2D.fa.gz out=read_misses2D.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow min=2 in=read_misses2D.fa.gz out=riboKmers2E.fa.gz fastawrap=4000 ow zl=9 pigz=16
bbduk.sh -Xmx8g in=read_misses2D.fa.gz ref=riboKmers2A.fa.gz,riboKmers2B.fa.gz,riboKmers2C.fa.gz,riboKmers2D.fa.gz,riboKmers2E.fa.gz out=read_misses2E.fa.gz zl=6 k=31 mm=f ordered fastawrap=4000 ow
kcompress.sh -Xmx31g ow in=riboKmers2A.fa.gz,riboKmers2B.fa.gz,riboKmers2C.fa.gz,riboKmers2D.fa.gz,riboKmers2E.fa.gz out=riboKmers2merged.fa.gz fastawrap=4000 ow zl=9 pigz=16
clumpify.sh k=16 in=riboKmers2merged.fa.gz out=riboKmers2clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers2clumped.fa.gz out=riboKmers2fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1

clumpify.sh k=16 in=riboKmers1A.fa.gz out=riboKmers1clumped.fa.gz g=1 zl=9 fastawrap=4000 reorder rcomp ow
fuse.sh -Xmx1g ow in=riboKmers1A.fa.gz out=riboKmers1fused.fa.gz fastawrap=8000 ow zl=11 pigz=32 maxlen=4000 npad=1


