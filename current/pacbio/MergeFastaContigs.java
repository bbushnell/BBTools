package pacbio;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Timer;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.Read;
import structures.ListNum;
import tracker.ReadStats;

/**
 * Merges FASTA or FASTQ records into synthetic chromosomes separated by blocks
 * of N bases. FASTA input may come from one file or a comma-separated list of
 * files; FASTQ input is converted to FASTA and treated one read per contig.
 * Records shorter than {@link #MIN_CONTIG_TO_ADD} are skipped, and output
 * chromosomes roll over when a pending contig plus padding would reach
 * {@link #MAX_OUTPUT_LEN}.
 * @author Brian Bushnell
 * @date Jul 10, 2012
 */
public class MergeFastaContigs{
	
	/**
	 * Parses merge settings, dispatches to the FASTA or FASTQ merge path, and
	 * prints global counters to standard output.
	 * Legacy positional fallbacks for {@code in}, {@code out}, {@code index},
	 * {@code npad}, and {@code npad2} are retained below, but ordinary
	 * positional arguments are rejected as unknown before the fallback runs.
	 * @param args Command-line arguments
	 */
	public static void main(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		String infile=null;
		String outfile=null;
		String outindex=null;
		int npl=-1;
		int npl2=-1;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(Parser.parseZip(arg, a, b)){
				//do nothing
			}else if(a.equals("in") && split.length>0){
				infile=b;
			}else if(a.equals("out") && split.length>0){
				outfile=b;
			}else if(a.equals("index") && split.length>0){
				outindex=b;
			}else if(a.equals("npad")){
				npl=N_PAD_LENGTH=Integer.parseInt(b);
			}else if(a.equals("npad2")){
				npl2=N_PAD_LENGTH2=Integer.parseInt(b);
			}else if(a.equals("maxdataout")){
				maxDataOut=Integer.parseInt(b);
			}else if(a.equals("mincontig")){
				MIN_CONTIG_TO_ADD=Integer.parseInt(b);
			}else if(a.equals("maxlen")){
				MAX_OUTPUT_LEN=Integer.parseInt(b);
			}else if(a.equals("maxchroms")){
				maxChromsOut=Integer.parseInt(b);
			}else if(a.equals("maxdata")){
				maxDataOut=Parse.parseKMG(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("append") || a.equals("app")){
				append=ReadStats.append=Parse.parseBoolean(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				overwrite=Parse.parseBoolean(b);
			}else if(a.equals("padfront") || a.equals("padstart")){
				PAD_START=Parse.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown argument "+arg);
			}
		}
		
		if(infile==null){infile=args[0];}
		if(outfile==null){outfile=args[1];}
		if(outindex==null){outindex=args[2];}
		
		try{
			if(npl<0 && args.length>3){N_PAD_LENGTH=Integer.parseInt(args[3]);}
			if(npl2<0 && args.length>4){N_PAD_LENGTH2=Integer.parseInt(args[4]);}
		}catch(NumberFormatException e){
			//ignore
		}
		
		if(infile.contains(".fq.") || infile.endsWith(".fq") || infile.contains(".fastq.") || infile.endsWith(".fastq")){
			mergeFastq(infile, outfile, outindex);
		}else{
			if(new File(infile).exists()){
//				System.err.println("Warning: This will run correctly, but I suggest against putting commas in your filenames.");
//				assert false : infile+", "+outfile+", "+outindex;
				mergeFasta(new String[] {infile}, outfile, outindex);
			}else{
				String[] files=infile.split(",");
				for(String s : files){
					if(!new File(s).exists()){throw new RuntimeException("Cannot find file "+s);}
				}
				mergeFasta(files, outfile, outindex);
			}
		}
		t.stop();
		
		System.out.println("MergeFastaContigs output for "+Arrays.toString(args));
		System.out.println("definedBasesIn:     \t"+definedBasesIn);
		System.out.println("contigsIn:          \t"+contigsIn);
		System.out.println("definedBasesOut:    \t"+definedBasesOut);
		System.out.println("basesOut:           \t"+dataOut);
		System.out.println("contigsOut:         \t"+contigsOut);
		System.out.println("chromsOut:          \t"+chromsOut);
		
		System.out.println("Time:\t"+t);
		
	}
	
	/**
	 * Legacy single-file FASTA merge path that writes merged FASTA to standard
	 * output and the contig index to standard error. The {@code outfile} and
	 * {@code outindex} parameters are retained for the historical signature but
	 * are not opened by this method; {@link #main(String[])} dispatches new FASTA
	 * work to {@link #mergeFasta(String[], String, String)}.
	 * @param infile Input FASTA file path
	 * @param outfile Ignored output FASTA path
	 * @param outindex Ignored output index path
	 */
	public static void merge(String infile, String outfile, String outindex){
		StringBuilder temp=new StringBuilder(MIN_CONTIG_TO_ADD);
		TextFile tf=new TextFile(infile, false);
		
//		OutputStream cos=ReadWrite.getOutputStream(outfile, false);
//		PrintWriter cpw=new PrintWriter(cos);
		
		long loc=N_PAD_LENGTH;
		int chrom=1;
		System.out.println(">chr"+chrom);
		npad=npad(N_PAD_LENGTH);
		printAsLines(npad, 0);
		
		String s=null;
		String label=null;
		for(s=tf.nextLine(); chrom<maxChromsOut && dataOut<maxDataOut; s=tf.nextLine()){
			if(s==null || s.charAt(0)=='>'){
				
				if(s!=null){contigsIn++;}
				
				//evict current contig
				if(temp.length()>=MIN_CONTIG_TO_ADD){
					
					long newloc=loc+temp.length()+N_PAD_LENGTH;
					if(newloc>=MAX_OUTPUT_LEN){
						//Evict old chrom

						//Make new chrom
						chrom++;
						loc=N_PAD_LENGTH;
						newloc=loc+temp.length()+N_PAD_LENGTH;
						System.out.println("\n>chr"+chrom);
						printAsLines(npad, 0);
					}
					
					printAsLines(temp, (int)(loc%lineBreak));
					
					definedBasesOut+=temp.length();
					contigsOut++;
					
					printAsLines(npad, (int)((loc+temp.length())%lineBreak));
					System.err.println(chrom+"\t"+loc+"\t"+label);
					loc=newloc;
				}else{
//					System.err.println("Ignored "+temp);
				}
				
				if(s==null){break;}
				temp.setLength(0);
				label=s.substring(1);
			}else{
				//append line to current contig
				temp.append(s);
				definedBasesIn+=s.length();
			}
		}
		tf.close();
		
		chromsOut=chrom;
		System.out.println();
	}
	
	/**
	 * Merges FASTA contigs from one or more files into synthetic chromosomes.
	 * The index records {@code chrom}, zero-based contig start, and original header
	 * for each emitted contig. Writers are started internally and joined through
	 * {@link TextStreamWriter#poisonAndWait()} after the final right pad is queued;
	 * the returned terminal-status booleans are ignored.
	 * @param infiles Input FASTA file paths
	 * @param outfile Output FASTA file path
	 * @param outindex Output index file path mapping contig starts
	 */
	public static void mergeFasta(String infiles[], String outfile, String outindex){
		
		if(new File(outfile).exists()){
			for(String s : infiles){assert(!s.equalsIgnoreCase(outfile));}
		}
		
		//if(verbose){System.err.println("A");}
		
		StringBuilder temp=new StringBuilder(MIN_CONTIG_TO_ADD);
		TextFile tf;

		TextStreamWriter cout=new TextStreamWriter(outfile, overwrite, false, false);
		TextStreamWriter iout=new TextStreamWriter(outindex, overwrite, false, false);
		
		cout.start();
		iout.start();
		//if(verbose){System.err.println("B");}
		
		long loc=(PAD_START ? N_PAD_LENGTH2 : 0);
		int chrom=1;
		cout.print(">chr"+chrom+"\n");
		npad=npad(N_PAD_LENGTH);
		npad2=npad2(N_PAD_LENGTH2);
		assert(npad.length()<=npad2.length());
		if(PAD_START){printAsLines(npad2, 0, cout);}
		boolean np2=true;
//		cout.poison();
//		assert(false) : "\n"+npad+"\n\n\n"+npad2+"\n";
		//if(verbose){System.err.println("C");}
		
//		assert(false) : PAD_START+", "+np2;
		
		for(String fname : infiles){
			tf=new TextFile(fname, false);
			String s=null;
			String label=null;
			if(verbose){System.err.println("Processing file "+fname);}
			for(s=tf.nextLine(); chrom<maxChromsOut && dataOut<maxDataOut; s=tf.nextLine()){
				//if(verbose){System.err.print("");}
				if(verbose){System.err.println("Processing line "+s);}
				if(s==null || s.charAt(0)=='>'){
					if(verbose){System.err.println("Contig break");}
//					System.err.println("chrom="+chrom+", maxChromsOut="+maxChromsOut);

					if(s!=null){contigsIn++;}
					if(verbose){System.err.println("Contigs="+contigsIn);}

					//evict current contig
					if(temp.length()>=MIN_CONTIG_TO_ADD){
						if(verbose){System.err.println("Big enough to add");}

						long newloc=loc+temp.length()+N_PAD_LENGTH;
						if(newloc>=MAX_OUTPUT_LEN){
							if(verbose){System.err.println("newloc>=MAX_OUTPUT_LEN");}
							//Evict old chrom
							printAsLines(npad2, (int)(loc%lineBreak), cout);

							//Make new chrom
							chrom++;
							loc=N_PAD_LENGTH2;
							newloc=loc+temp.length()+N_PAD_LENGTH;
							cout.print("\n>chr"+chrom+"\n");
							if(PAD_START){printAsLines(npad2, 0, cout);}
							np2=true;
						}
						if(verbose){System.err.println("G");}

						printAsLines(temp, (int)(loc%lineBreak), cout);

						definedBasesOut+=temp.length();
						contigsOut++;
						
						if(np2){
							if(verbose){System.err.println("np2");}
							if(PAD_START){
								if(verbose){System.err.println("PAD_START");}
								loc=N_PAD_LENGTH2;
								newloc=N_PAD_LENGTH2+temp.length();
							}else{
								if(verbose){System.err.println("~PAD_START");}
								loc=0;
								newloc=temp.length();
							}
						}else{
							if(verbose){System.err.println("PAD_START");}
							printAsLines(npad, (int)((loc+temp.length())%lineBreak), cout);
						}
						if(verbose){System.err.println("H");}
						if(label!=null){iout.print(chrom+"\t"+loc+"\t"+label+"\n");}
						loc=newloc;
						np2=false;
					}else{
						//					System.err.println("Ignored "+temp);
					}
					if(verbose){System.err.println("Done with contig");}
					
					temp.setLength(0);
					if(s==null){break;}
					label=s.substring(1);
				}else{
					np2=false;
					//if(verbose){System.err.print("J");}
					//append line to current contig
					temp.append(s);
					definedBasesIn+=s.length();
					if(verbose){System.err.println("Normal line.  definedBasesIn="+definedBasesIn);}
				}
				//if(verbose){System.err.print("K");}
			}
			tf.close();
			//if(verbose){System.err.print("L");}
		}
		//if(verbose){System.err.println("M");}
		
		chromsOut=chrom;
		
		assert(temp.length()==0) : temp.length();
		printAsLines(npad2, (int)(loc%lineBreak), cout);
		//if(verbose){System.err.println("N");}
		
		
		cout.print("\n");
		cout.poisonAndWait();
		iout.poisonAndWait();
		//if(verbose){System.err.println("O");}
	}
	
	
	
	/**
	 * Converts a FASTQ file to padded synthetic FASTA chromosomes.
	 * Each read is treated as one candidate contig and recorded in the index by
	 * read ID. The input stream is closed explicitly; output writers are poisoned
	 * without waiting for completion or checking their terminal status.
	 * @param in1 Input FASTQ file path
	 * @param outfile Output FASTA file path
	 * @param outindex Output index file path mapping read starts
	 */
	public static void mergeFastq(String in1, String outfile, String outindex){
		StringBuilder temp=new StringBuilder(MIN_CONTIG_TO_ADD);
		
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.DETECT_QUALITY=false;
		long maxReads=-1;
		
		final ConcurrentReadInputStream cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ff1, null);
//			if(verbose){System.err.println("Started cris");}
			cris.start(); //4567
		}
		

		TextStreamWriter cout=new TextStreamWriter(outfile, overwrite, false, false);
		TextStreamWriter iout=new TextStreamWriter(outindex, overwrite, false, false);
		
		cout.start();
		iout.start();
		
		long loc=N_PAD_LENGTH2;
		int chrom=1;
		cout.print(">chr"+chrom+"\n");
		npad=npad(N_PAD_LENGTH);
		npad2=npad2(N_PAD_LENGTH2);
		assert(npad.length()<=npad2.length());
		printAsLines(npad2, 0, cout);
		
		
		String s=null;
		String label=null;
		
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);


		while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
			//System.err.println("reads.size()="+reads.size());
			for(Read r : reads){

				s=new String(r.bases);
				label=r.id;

				temp.append(s);

				if(temp.length()>=MIN_CONTIG_TO_ADD){

					long newloc=loc+temp.length()+N_PAD_LENGTH;
					if(newloc>=MAX_OUTPUT_LEN){
						//Evict old chrom
						printAsLines(npad2, (int)(loc%lineBreak), cout);

						//Make new chrom
						chrom++;
						loc=N_PAD_LENGTH2;
						newloc=loc+temp.length()+N_PAD_LENGTH;
						cout.print("\n>chr"+chrom+"\n");
						printAsLines(npad2, 0, cout);
					}

					printAsLines(temp, (int)(loc%lineBreak), cout);
					printAsLines(npad, (int)((loc+temp.length())%lineBreak), cout);
					iout.println(chrom+"\t"+loc+"\t"+label);
					loc=newloc;
				}else{
					//						System.err.println("Ignored "+temp);
				}

				temp.setLength(0);
				if(s==null){break;}
				label=s.substring(1);


			}
			//System.err.println("returning list");
			cris.returnList(ln);
			//System.err.println("fetching list");
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}

		
		cris.returnList(ln);
		
		assert(temp.length()==0) : temp.length();
		printAsLines(npad2, (int)(loc%lineBreak), cout);
		

		ReadWrite.closeStream(cris);
		
		cout.print("\n");
		cout.poison();
		iout.poison();
	}
	
	/**
	 * Prints sequence to standard output with line breaks aligned to lineBreak.
	 * TODO: Probable legacy stdout bug: {@code mod==0} starts after the first
	 * chunk, unlike {@link #printAsLines(CharSequence, int, TextStreamWriter)}.
	 * @param sb Sequence or padding to print
	 * @param mod Existing bases on the current output line
	 */
	private static void printAsLines(CharSequence sb, int mod){
		dataOut+=sb.length();
		assert(mod<lineBreak);
		if(mod>0){
			CharSequence s=sb.subSequence(0, min(lineBreak-mod, sb.length()));
			if(s.length()+mod==lineBreak){
				System.out.println(s);
			}else{
				System.out.print(s);
			}
		}
		
		int loc=lineBreak-mod;
		for(; loc<sb.length(); loc+=lineBreak){
			CharSequence s=sb.subSequence(loc, min(loc+lineBreak, sb.length()));
			if(s.length()==lineBreak){
				System.out.println(s);
			}else{
				System.out.print(s);
			}
		}
	}

	/**
	 * Queues sequence to a writer with line breaks aligned to lineBreak.
	 * @param sb Sequence or padding to print
	 * @param mod Existing bases on the current output line
	 * @param cout Destination FASTA writer
	 */
	private static void printAsLines(CharSequence sb, int mod, TextStreamWriter cout){
		dataOut+=sb.length();
		assert(mod<lineBreak);
		if(mod>0){
			
			CharSequence s=sb.subSequence(0, min(lineBreak-mod, sb.length()));
			
//			System.out.println(mod+", "+s.length()+", "+(s.length()+mod)+", "+lineBreak);
			
			if(s.length()+mod==lineBreak){
				cout.println(s);
			}else{
				cout.print(s);
			}
		}
		
		int loc=(mod==0 ? 0 : lineBreak-mod);
		for(; loc<sb.length(); loc+=lineBreak){
			CharSequence s=sb.subSequence(loc, min(loc+lineBreak, sb.length()));
			
//			System.out.println(mod+", "+s.length()+", "+(s.length()+mod)+", "+lineBreak+", loc="+loc);
			
			if(s.length()==lineBreak){
				cout.println(s);
			}else{
				cout.print(s);
			}
		}
	}

	/**
	 * Returns and caches an N block of the requested inter-contig pad length.
	 * @param N_PAD_LENGTH Number of N bases
	 * @return Cached N block
	 */
	public static String npad(int N_PAD_LENGTH){
		if(npad==null || npad.length()!=N_PAD_LENGTH){
			StringBuilder sb=new StringBuilder(N_PAD_LENGTH);
			for(int i=0; i<N_PAD_LENGTH; i++){
				sb.append('N');
			}
			npad=sb.toString();
		}
		return npad;
	}

	/**
	 * Returns and caches an N block of the requested end-pad length.
	 * @param N_PAD_LENGTH Number of N bases
	 * @return Cached N block
	 */
	public static String npad2(int N_PAD_LENGTH){
		if(npad2==null || npad2.length()!=N_PAD_LENGTH){
			StringBuilder sb=new StringBuilder(N_PAD_LENGTH);
			for(int i=0; i<N_PAD_LENGTH; i++){
				sb.append('N');
			}
			npad2=sb.toString();
		}
		return npad2;
	}
	
	private static final int min(int x, int y){return x<y ? x : y;}
	@SuppressWarnings("unused")
	private static final int max(int x, int y){return x>y ? x : y;}
	
	static long definedBasesIn=0;
	static long contigsIn=0;
	static long definedBasesOut=0;
	static long contigsOut=0;
	static long chromsOut=0;
	
	/** Bases per FASTA output line */
	public static int lineBreak=80;
	/** Legacy unused lineBreak-derived field */
	public static int modulo=lineBreak+1;
	/** N bases between emitted contigs */
	public static int N_PAD_LENGTH=300;
	/** N bases used for chromosome-end padding */
	public static int N_PAD_LENGTH2=2000;
	/** Shorter contigs are skipped */
	public static int MIN_CONTIG_TO_ADD=150;
	/**
	 * Maximum requested chromosome length before rollover; the 220M default
	 * preserves a historical capacity target.
	 */
	public static int MAX_OUTPUT_LEN=220000000;
	/** Maximum number of synthetic chromosomes */
	public static int maxChromsOut=60000;
	/** Approximate maximum number of queued sequence and padding bases */
	public static long maxDataOut=Long.MAX_VALUE;
	/** Count of queued sequence and padding bases, excluding headers and newlines */
	private static long dataOut=0;
	/** Cached N-padding strings */
	public static String npad, npad2;
	/** Overwrite existing output files */
	public static boolean overwrite=true;
	/** Append to existing output files */
	public static boolean append=false;
	/** Add end padding before the first contig */
	public static boolean PAD_START=true;
	/** Print per-line FASTA progress */
	public static boolean verbose=false;
	
}
