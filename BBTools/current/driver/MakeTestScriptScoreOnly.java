package driver;

public class MakeTestScriptScoreOnly {
	
	
	public static void main(String[] args){
		
		assert(args.length==1) : "Please enter number of reads.";
		numReads=Integer.parseInt(args[0]);
		
//		String[] strings=new String[] {
//				"/work/bbushnell/ssaha2/ssaha2 -solexa -outfile #S.sam -best -1 -output sam_soft -save hg37 " +
//				"/work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime defaultTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S.sam #R ssaha2",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime defaultTime.txt"
//		};
		
		
//		String[] strings=new String[] {
//				"bwa aln -t 22 bs_ /work/bbushnell/synth/reads_B1_100000x100bp_#S_chr1-25.fq > temp_default.sai",
//				"bwa samse bs_ temp_default.sai /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq > #S_default.sam",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime defaultTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_default.sam #R",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime defaultTime.txt"
//		};
		
		
//		String[] strings=new String[] {
//				"java -ea -Xms24g -Xmx31g -server -XX:+UseNUMA -XX:+AggressiveOpts -XX:+UseCompressedOops " +
//				"align.TestIndex11f 1 25 100 0 /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq null " +
//				"outfile=#S_bbmap11f.sam cs=false threads=22 paired=false pairlen=100 build=37 match=short " +
//				"removeambiguous=false fastqparsecustom overwrite savepar=false",
//
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bbmap11fTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_bbmap11f.sam #R",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bbmap11fTime.txt"
//		};
		
		
//		String[] strings=new String[] {
//				"bowtie --best -y --chunkmbs 1024 --strata -m 1 -k 2 -v 3 -p 24 -t -q -S HG37" +
//				" /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq #S_bowtie.sam",
//
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bowtieTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_bowtie.sam #R",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bowtieTime.txt"
//		};
		
		
//		String[] strings=new String[] {
//				"bfast match -T $TMPDIR/ -n 16 -f hg19.fa -r /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq > $TMPDIR/#S.bmf",
//				"bfast localalign -n 16 -f hg19.fa -m $TMPDIR/#S.bmf > $TMPDIR/#S.baf",
////				"bfast postprocess -n 16 -a 3 -f hg19.fa -i $TMPDIR/#S.baf > #S.sam",
////				"bfast postprocess -n 16 -a 3 -m 20 -f hg19.fa -i $TMPDIR/#S.baf > #S_r#R.sam",
//				"bfast postprocess -n 16 -M 20 -f hg19.fa -i $TMPDIR/#S.baf > #S_r#R.sam",
//
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bfastTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_r#R.sam #R",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime bfastTime.txt"
//		};
		
//		String[] strings=new String[] {
//				"smalt_x86_64 map -n 8 -a -f samsoft -o #S_r#R.sam hg37 /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_r#R.sam #R ssaha2",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
//		};
		
//		String[] strings=new String[] {
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
//				"./soap -p 24 -a /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq -D hg37.fa.index -o #S_r#R.soap",
//				"perl soap2sam.pl -p #S_r#R.soap > #S_r#R.sam",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
//				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_r#R.sam #R",
//		};
		
		String[] strings=new String[] {
				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
				"./bin/gmapper-ls /work/bbushnell/synth/reads_B1_#Rx100bp_#S_chr1-25.fq --single-best-mapping --qv-offset 33 -L hg37 -N 24 -o 5 -h 80% > #S_r#R.sam",
				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.PrintTime mapTime.txt",
				"java -ea -Xmx96m -cp /work/bbushnell/java/ align.GradeSamFile #S_r#R.sam #R",
		};
		
		int[] blank=new int[] {0, 0, 0, 0};

		print(strings, blank, 100);
		print(strings, blank, 100);
		print(strings, blank, 100);
		print(strings, blank, 100);
		for(int[] array : sets){
			print(strings, array, numReads);
		}
		
	}
	
	private static void print(String[] array, int[] blank, int x) {
		
		String counts=(blank[0]+"S_"+blank[1]+"I_"+blank[2]+"D_"+blank[3]+"U");
		String reads=""+x;
		
		for(String s : array){
			String s2=s.replaceAll("#S", counts).replaceAll("#R", reads);
			System.out.println(s2);
		}
		System.out.println();
		
	}

	public static int numReads=400000;
	
	public static final int[][] sets=new int[][] {
		{0, 0, 0, 0},
		{1, 0, 0, 0},
		{2, 0, 0, 0},
		{3, 0, 0, 0},
		{4, 0, 0, 0},
		{5, 0, 0, 0},
		{6, 0, 0, 0},
		{7, 0, 0, 0},
		{8, 0, 0, 0},
		{10, 0, 0, 0},
		{12, 0, 0, 0},
		{14, 0, 0, 0},
		{16, 0, 0, 0},
		{18, 0, 0, 0},
		{20, 0, 0, 0},
		{24, 0, 0, 0},
		{28, 0, 0, 0},
		{32, 0, 0, 0},
		{36, 0, 0, 0},
		{40, 0, 0, 0},

		{0, 1, 0, 0},
		{0, 2, 0, 0},
		{0, 3, 0, 0},
		{0, 4, 0, 0},
		{0, 5, 0, 0},
		{0, 6, 0, 0},
		{0, 7, 0, 0},
		{0, 8, 0, 0},
		{0, 10, 0, 0},
		{0, 12, 0, 0},
		{0, 14, 0, 0},
		{0, 16, 0, 0},
		{0, 20, 0, 0},
		{0, 24, 0, 0},
		{0, 28, 0, 0},
		{0, 32, 0, 0},
		{0, 36, 0, 0},
		{0, 40, 0, 0},

		{0, 0, 1, 0},
		{0, 0, 2, 0},
		{0, 0, 3, 0},
		{0, 0, 4, 0},
		{0, 0, 5, 0},
		{0, 0, 6, 0},
		{0, 0, 7, 0},
		{0, 0, 8, 0},
		{0, 0, 10, 0},
		{0, 0, 12, 0},
		{0, 0, 14, 0},
		{0, 0, 16, 0},
		{0, 0, 20, 0},
		{0, 0, 24, 0},
		{0, 0, 28, 0},
		{0, 0, 32, 0},
		{0, 0, 48, 0},
		{0, 0, 64, 0},
		{0, 0, 128, 0},
		{0, 0, 192, 0},
		{0, 0, 256, 0},
		{0, 0, 512, 0},
		{0, 0, 1000, 0},
		{0, 0, 2000, 0},
		{0, 0, 3000, 0},
		{0, 0, 4000, 0},
		{0, 0, 6000, 0},
		{0, 0, 8000, 0},
		{0, 0, 10000, 0},
		{0, 0, 12000, 0},
		{0, 0, 14000, 0},
		{0, 0, 16000, 0},
		{0, 0, 20000, 0},
		{0, 0, 24000, 0},
		{0, 0, 28000, 0},
		{0, 0, 32000, 0},

		{0, 0, 0, 1},
		{0, 0, 0, 2},
		{0, 0, 0, 3},
		{0, 0, 0, 4},
		{0, 0, 0, 5},
		{0, 0, 0, 6},
		{0, 0, 0, 7},
		{0, 0, 0, 8},
		{0, 0, 0, 10},
		{0, 0, 0, 12},
		{0, 0, 0, 14},
		{0, 0, 0, 16},
		{0, 0, 0, 20},
		{0, 0, 0, 24},
		{0, 0, 0, 28},
		{0, 0, 0, 32},
		{0, 0, 0, 36},
		{0, 0, 0, 40}
	};
	
}
