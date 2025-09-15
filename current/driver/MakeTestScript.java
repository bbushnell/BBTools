package driver;

public class MakeTestScript {
	
	
	public static void main(String[] args){
		
		assert(args.length>=1) : "Please enter number of reads.";
		numReads=Integer.parseInt(args[0]);
		readlen=Integer.parseInt(args[1]);
		
		String mode=args[2];
		String extra=(args.length>3 ? args[3] : "EXTRA");
		
		String printtime="java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime ";
		String gradesam="java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.GradeSamFile ";
		String time=mode+"Time.txt";
		
		String[] strings=null;
		
//		strings=new String[] {
//				"/house/homedirs/b/bushnell/ssaha2/ssaha2 -solexa -outfile #S.sam -best -1 -output sam_soft -save hg37 " +
//				"reads_B1_#Rx#Lbp_#S.fastq",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime defaultTime.txt",
//				gradesam+"#S.sam #R ssaha2",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime defaultTime.txt"
//		};
		
		if(mode.equalsIgnoreCase("bwa")){
			strings=new String[] {
//				printtime+time+" false",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa aln -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fastq 1>temp_bwa.sai",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa samse "+extra+" temp_bwa.sai reads_B1_#Rx#Lbp_#S.fastq 1>bwa_#S_r#Rx#L.sam",
//				printtime+time,
//				gradesam+"bwa_#S_r#Rx#L.sam #R",
				
				printtime+time+" false",
				"/house/homedirs/b/bushnell/bwa/bwa aln -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fastq 1>bwa_#S_r#Rx#L.sai",
				"/house/homedirs/b/bushnell/bwa/bwa samse "+extra+" bwa_#S_r#Rx#L.sai reads_B1_#Rx#Lbp_#S.fastq 1>bwa_#S_r#Rx#L.sam",
				printtime+time,
				gradesam+"bwa_#S_r#Rx#L.sam #R",
			};
		}
		
		if(mode.equalsIgnoreCase("bwamem")){
			strings=new String[] {
//				printtime+time+" false",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa aln -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fastq 1>temp_bwa.sai",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa samse "+extra+" temp_bwa.sai reads_B1_#Rx#Lbp_#S.fastq 1>bwa_#S_r#Rx#L.sam",
//				printtime+time,
//				gradesam+"bwa_#S_r#Rx#L.sam #R",
				
				printtime+time+" false",
				"/house/homedirs/b/bushnell/bwa74/bwa mem -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fastq 1>bwamem_#S_r#Rx#L.sam",
				printtime+time,
				gradesam+"bwamem_#S_r#Rx#L.sam #R",
			};
		}
		
		if(mode.equalsIgnoreCase("bwasw")){
			strings=new String[] {
//				printtime+time+" false",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa aln -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fastq 1>temp_bwa.sai",
//				"memtime /house/homedirs/b/bushnell/bwa/bwa samse "+extra+" temp_bwa.sai reads_B1_#Rx#Lbp_#S.fastq 1>bwa_#S_r#Rx#L.sam",
//				printtime+time,
//				gradesam+"bwa_#S_r#Rx#L.sam #R",
				
				printtime+time+" false",
				"/house/homedirs/b/bushnell/bwa/bwa bwasw -b5 -q2 -r1 -z10 -t 32 "+extra+" reads_B1_#Rx#Lbp_#S.fasta 1>bwa_#S_r#Rx#L.sam",
				printtime+time,
				gradesam+"bwa_#S_r#Rx#L.sam #R",
			};
		}
		
		if(mode.startsWith("bbmap")){
			int k=13;
			String s2=mode.replaceFirst("bbmap", "");
			if(s2.length()>0){
				k=Integer.parseInt(s2);
			}
			strings=new String[] {
					printtime+time+" false",
				"memtime java -ea -Xmx106g -cp /house/homedirs/b/bushnell/beta18/ " +
				"align2.BBMap in=reads_B1_#Rx#Lbp_#S.fastq out=bbmap"+k+"_#S_r#Rx#L.sam overwrite k="+k+" printtoerr",
				printtime+time,
				gradesam+"bbmap"+k+"_#S_r#Rx#L.sam #R",
			};
		}
		
		if(mode.equalsIgnoreCase("bowtie2")){
			strings=new String[] {
					printtime+time+" false",
				"memtime bowtie2 -x bow2ref -U reads_B1_#Rx#Lbp_#S.fastq -S bowtie2_#S_r#Rx#L.sam --phred33 -p 32",
				printtime+time,
				gradesam+"bowtie2_#S_r#Rx#L.sam #R",
			};
		}
		
		if(mode.equalsIgnoreCase("gsnap")){
			strings=new String[] {
					printtime+time+" false",
				"memtime /house/homedirs/b/bushnell/gsnap/bin/gsnap -t 32 -d "+extra+" -A sam reads_B1_#Rx#Lbp_#S.fastq > gsnap_#S_r#Rx#L.sam",
				printtime+time,
				gradesam+"gsnap_#S_r#Rx#L.sam #R",
			};
		}
		
		
//		strings=new String[] {
//				"bowtie --best -y --chunkmbs 1024 --strata -m 1 -k 2 -v 3 -p 24 -t -q -S HG37" +
//				" reads_B1_#Rx#Lbp_#S.fastq #S_bowtie.sam",
//
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime bowtieTime.txt",
//				gradesam+"#S_bowtie.sam #R",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime bowtieTime.txt"
//		};
		
		
//		strings=new String[] {
//				"bfast match -T $TMPDIR/ -n 16 -f hg19.fa -r reads_B1_#Rx#Lbp_#S.fastq > $TMPDIR/#S.bmf",
//				"bfast localalign -n 16 -f hg19.fa -m $TMPDIR/#S.bmf > $TMPDIR/#S.baf",
////				"bfast postprocess -n 16 -a 3 -f hg19.fa -i $TMPDIR/#S.baf > #S.sam",
////				"bfast postprocess -n 16 -a 3 -m 20 -f hg19.fa -i $TMPDIR/#S.baf > #S_r#Rx#L.sam",
//				"bfast postprocess -n 16 -M 20 -f hg19.fa -i $TMPDIR/#S.baf > #S_r#Rx#L.sam",
//
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime bfastTime.txt",
//				gradesam+"#S_r#Rx#L.sam #R",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime bfastTime.txt"
//		};

		if(mode.equalsIgnoreCase("smalt")){
			strings=new String[] {
					printtime+time+" false",
				"memtime /house/homedirs/b/bushnell/smalt/smalt_x86_64 map -n 32 -f sam -o smalt_#S_r#Rx#L.sam smaltindex reads_B1_#Rx#Lbp_#S.fastq",
				printtime+time,
				gradesam+"smalt_#S_r#Rx#L.sam #R ssaha2",
			};
		}

		if(mode.equalsIgnoreCase("snap")){
			strings=new String[] {
					printtime+time+" false",
				"memtime /house/homedirs/b/bushnell/snap/snap single snapref reads_B1_#Rx#Lbp_#S.fastq -o snap_#S_r#Rx#L.sam -t 32 -b",
				printtime+time,
				gradesam+"snap_#S_r#Rx#L.sam #R",
			};
		}

		if(mode.equalsIgnoreCase("masai")){
			strings=new String[] {
					printtime+time+" false",
				"memtime /house/homedirs/b/bushnell/masai/masai_mapper --output-format sam "+extra+" reads_B1_#Rx#Lbp_#S.fastq",
				printtime+time,
				gradesam+"reads_B1_#Rx#Lbp_#S.sam #R",
			};
		}

		if(mode.equalsIgnoreCase("blasr")){
			System.out.println("source /house/sdm/pacbio/smrtanalysis-installs/smrtanalysis-2.0.0/etc/setup.sh\n");
			strings=new String[] {
					printtime+time+" false",
				"memtime blasr reads_B1_#Rx#Lbp_#S.fastq "+extra+" -sam -out blasr_#S_r#Rx#L.sam -bestn 1 -nproc 32",
				printtime+time,
				gradesam+"blasr_#S_r#Rx#L.sam #R blasr",
			};
		}
		
//		strings=new String[] {
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				"./soap -p 24 -a reads_B1_#Rx#Lbp_#S.fastq -D hg37.fa.index -o #S_r#Rx#L.soap",
//				"perl soap2sam.pl -p #S_r#Rx#L.soap > #S_r#Rx#L.sam",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				gradesam+"#S_r#Rx#L.sam #R",
//		};
		
//		strings=new String[] {
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				"./bin/gmapper-ls reads_B1_#Rx#Lbp_#S.fastq --single-best-mapping --qv-offset 33 -L hg37 -N 24 -o 5 -h 80% > #S_r#Rx#L.sam",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				gradesam+"#S_r#Rx#L.sam #R",
//		};
		
//		strings=new String[] {
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				"./bin/MosaikBuild -q reads_B1_#Rx#Lbp_#S.fastq -out $TMPDIR/reads_B1_#Rx100bp_#S_chr1-25.dat -st illumina",
//				"./bin/MosaikAligner -in $TMPDIR/reads_B1_#Rx100bp_#S_chr1-25.dat -out $TMPDIR/reads_B1_#Rx100bp_#S_chr1-25_aligned.dat -ia hg37_ref.dat -hs 15 -bw=29 -j hg37_jumpdb -act 20 -mm 32 -mhp 100 -p 32 -m unique",
//				"./bin/MosaikText -in $TMPDIR/reads_B1_#Rx100bp_#S_chr1-25_aligned.dat -sam #S_r#Rx#L.sam",
//				"java -ea -Xmx96m -cp /house/homedirs/b/bushnell/beta18/ align2.PrintTime mapTime.txt",
//				gradesam+"#S_r#Rx#L.sam #R",
//		};
		
		int[] blank=new int[] {0, 0, 0, 0, 0};
		
		int preload=100;
		if(mode.equalsIgnoreCase("masai")){
			preload=1000;
		}
		print(strings, blank, preload);
		print(strings, blank, preload);
		print(strings, blank, preload);
		print(strings, blank, preload);
		for(int[] array : sets){
			print(strings, array, numReads);
		}
		
	}
	
	private static void print(String[] array, int[] blank, int x) {
		
		int rl=readlen;
		if(blank.length>5){rl=blank[5];}
		
		String counts=(blank[0]+"S_"+blank[1]+"I_"+blank[2]+"D_"+blank[3]+"U_"+blank[4]+"N");
		String reads=""+x;
		String len=""+rl;
		
		for(String s : array){
			String s2=s.replaceAll("#S", counts).replaceAll("#R", reads).replaceAll("#L", len);
			System.out.println(s2);
		}
		System.out.println();
		
	}

	public static int numReads=400000;
	public static int readlen=150;
	
	public static final int[][] sets=new int[][] {
		{0, 0, 0, 0, 0},
		{1, 0, 0, 0, 0},
		{2, 0, 0, 0, 0},
		{3, 0, 0, 0, 0},
		{4, 0, 0, 0, 0},
		{5, 0, 0, 0, 0},
		{6, 0, 0, 0, 0},
		{7, 0, 0, 0, 0},
		{8, 0, 0, 0, 0},
		{10, 0, 0, 0, 0},
		{12, 0, 0, 0, 0},
		{14, 0, 0, 0, 0},
		{16, 0, 0, 0, 0},
		{18, 0, 0, 0, 0},
		{20, 0, 0, 0, 0},
		{24, 0, 0, 0, 0},
		{28, 0, 0, 0, 0},
		{32, 0, 0, 0, 0},
		{36, 0, 0, 0, 0},
		{40, 0, 0, 0, 0},

		{0, 1, 0, 0, 0},
		{0, 2, 0, 0, 0},
		{0, 3, 0, 0, 0},
		{0, 4, 0, 0, 0},
		{0, 5, 0, 0, 0},
		{0, 6, 0, 0, 0},
		{0, 7, 0, 0, 0},
		{0, 8, 0, 0, 0},
		{0, 10, 0, 0, 0},
		{0, 12, 0, 0, 0},
		{0, 14, 0, 0, 0},
		{0, 16, 0, 0, 0},
		{0, 18, 0, 0, 0},
		{0, 20, 0, 0, 0},
		{0, 24, 0, 0, 0},
		{0, 28, 0, 0, 0},
		{0, 32, 0, 0, 0},
		{0, 36, 0, 0, 0},
		{0, 40, 0, 0, 0},

		{0, 0, 1, 0, 0},
		{0, 0, 2, 0, 0},
		{0, 0, 3, 0, 0},
		{0, 0, 4, 0, 0},
		{0, 0, 5, 0, 0},
		{0, 0, 6, 0, 0},
		{0, 0, 7, 0, 0},
		{0, 0, 8, 0, 0},
		{0, 0, 10, 0, 0},
		{0, 0, 12, 0, 0},
		{0, 0, 14, 0, 0},
		{0, 0, 16, 0, 0},
		{0, 0, 18, 0, 0},
		{0, 0, 20, 0, 0},
		{0, 0, 24, 0, 0},
		{0, 0, 28, 0, 0},
		{0, 0, 32, 0, 0},
		{0, 0, 36, 0, 0},
		{0, 0, 40, 0, 0},
		{0, 0, 48, 0, 0},
		{0, 0, 56, 0, 0},
		{0, 0, 64, 0, 0},
		{0, 0, 96, 0, 0},
		{0, 0, 128, 0, 0},
		{0, 0, 192, 0, 0},
		{0, 0, 256, 0, 0},
		{0, 0, 384, 0, 0},
		{0, 0, 512, 0, 0},
		{0, 0, 768, 0, 0},
		{0, 0, 1000, 0, 0},
		{0, 0, 1500, 0, 0},
		{0, 0, 2000, 0, 0},
		{0, 0, 3000, 0, 0},
		{0, 0, 4000, 0, 0},
		{0, 0, 6000, 0, 0},
		{0, 0, 8000, 0, 0},
		{0, 0, 12000, 0, 0},
		{0, 0, 16000, 0, 0},
		{0, 0, 24000, 0, 0},
		{0, 0, 32000, 0, 0},
		{0, 0, 48000, 0, 0},
		{0, 0, 64000, 0, 0},
		{0, 0, 96000, 0, 0},
		{0, 0, 128000, 0, 0},

		{0, 0, 0, 1, 0},
		{0, 0, 0, 2, 0},
		{0, 0, 0, 3, 0},
		{0, 0, 0, 4, 0},
		{0, 0, 0, 5, 0},
		{0, 0, 0, 6, 0},
		{0, 0, 0, 7, 0},
		{0, 0, 0, 8, 0},
		{0, 0, 0, 10, 0},
		{0, 0, 0, 12, 0},
		{0, 0, 0, 14, 0},
		{0, 0, 0, 16, 0},
		{0, 0, 0, 18, 0},
		{0, 0, 0, 20, 0},
		{0, 0, 0, 24, 0},
		{0, 0, 0, 28, 0},
		{0, 0, 0, 32, 0},
		{0, 0, 0, 36, 0},
		{0, 0, 0, 40, 0},

		{0, 0, 0, 0, 1},
		{0, 0, 0, 0, 2},
		{0, 0, 0, 0, 3},
		{0, 0, 0, 0, 4},
		{0, 0, 0, 0, 5},
		{0, 0, 0, 0, 6},
		{0, 0, 0, 0, 7},
		{0, 0, 0, 0, 8},
		{0, 0, 0, 0, 10},
		{0, 0, 0, 0, 12},
		{0, 0, 0, 0, 14},
		{0, 0, 0, 0, 16},
		{0, 0, 0, 0, 18},
		{0, 0, 0, 0, 20},
		{0, 0, 0, 0, 24},
		{0, 0, 0, 0, 28},
		{0, 0, 0, 0, 32},
		{0, 0, 0, 0, 36},
		{0, 0, 0, 0, 40},

		{0, 0, 0, 0, 0, 400},
		{2, 2, 2, 2, 0, 400},
		{4, 2, 2, 2, 0, 400},
		{6, 3, 3, 3, 0, 400},
		{8, 4, 4, 4, 0, 400},
		{10, 4, 4, 4, 0, 400},
		{12, 4, 4, 4, 0, 400},
		{14, 4, 4, 4, 0, 400},
		{16, 4, 4, 4, 0, 400},
		{18, 4, 4, 4, 0, 400},
		{20, 5, 5, 5, 0, 400},
	};
	
}
