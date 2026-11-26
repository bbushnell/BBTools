package test.jgi;

import java.io.File;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Generates the deterministic CalcTrueQuality q12 timing input and independent
 * q10/old-q12/next-q12/p oracles by scaling the accepted JGI-178 three-record SAM.
 * This writes input data only; it does not run CalcTrueQuality or touch product classes.
 * @author Jean
 */
public final class CalcTrueQualityQ12TimingData{

	/** Generates one bounded triplet count. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CalcTrueQualityQ12TimingData.class, false);
		try{
			final Args parsed=new Args(pp.args);
			final CalcTrueQualityQ12TimingData maker=
				new CalcTrueQualityQ12TimingData(parsed.triplets, parsed.overwrite,
					parsed.sam, parsed.q10, parsed.q12old, parsed.q12next, parsed.p);
			maker.process();
			System.out.print("CALCTRUEQUALITY_Q12_TIMING_DATA_PASS triplets="+maker.triplets+
				" samBytes="+maker.samBytes+" q10Bytes="+maker.q10Bytes+
				" q12oldBytes="+maker.q12oldBytes+" q12nextBytes="+maker.q12nextBytes+
				" pBytes="+maker.pBytes+"\n");
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** Stores validated arguments and checks the source-level byte-bound table. */
	private CalcTrueQualityQ12TimingData(final int triplets_, final boolean overwrite_,
			final String sam_, final String q10_, final String q12old_,
			final String q12next_, final String p_){
		triplets=triplets_;
		overwrite=overwrite_;
		sam=sam_;
		q10=q10_;
		q12old=q12old_;
		q12next=q12next_;
		p=p_;
		assert(checkByteBounds()) : "The finite JGI-178 q12 timing byte-bound table "+
			"must match this generator's literal output.";
	}

	/** Writes the SAM and four oracle files, then verifies exact uncompressed file lengths. */
	private void process(){
		validate();
		samBytes=writeSam();
		q10Bytes=writeQualityOracle(q10, OLD_Q12_ROWS, expectedQualityBytes(OLD_Q12_ROWS, triplets));
		q12oldBytes=writeQualityOracle(q12old, OLD_Q12_ROWS, expectedQualityBytes(OLD_Q12_ROWS, triplets));
		q12nextBytes=writeQualityOracle(q12next, NEXT_Q12_ROWS, expectedQualityBytes(NEXT_Q12_ROWS, triplets));
		pBytes=writeCycleOracle(p, P_ROWS, expectedCycleBytes(P_ROWS, triplets));
	}

	/** Checks bounded count, output uniqueness, uncompressed paths and overwrite policy. */
	private void validate(){
		require(isAllowedTriplets(triplets), "triplets must be one of 1, 250k, 1m, 4m, 8m; observed "+triplets);
		final String[] names={sam, q10, q12old, q12next, p};
		final String[] labels={"out", "q10", "q12old", "q12next", "p"};
		for(int i=0; i<names.length; i++){
			require(names[i]!=null && names[i].length()>0, "Missing required "+labels[i]+"=<path>.");
			require(!names[i].equalsIgnoreCase("stdout") && !names[i].startsWith("/dev/std"),
				labels[i]+" must name a regular uncompressed file: "+names[i]);
			require(!ReadWrite.isCompressed(names[i]), labels[i]+" must be uncompressed for exact byte bounds: "+names[i]);
			final File file=new File(names[i]);
			final File parent=file.getAbsoluteFile().getParentFile();
			require(parent!=null && parent.isDirectory(), labels[i]+" parent directory does not exist: "+names[i]);
			require(!file.isDirectory(), labels[i]+" points to a directory: "+names[i]);
			require(overwrite || !file.exists(), labels[i]+" already exists and ow=f: "+names[i]);
			for(int j=0; j<i; j++){
				require(!sameFile(file, new File(names[j])), labels[i]+" reuses "+labels[j]+": "+names[i]);
			}
		}
	}

	/** Writes the scaled SAM with unique qnames and accepted non-qname fields. */
	private long writeSam(){
		final long expected=expectedSamBytes(triplets);
		ByteStreamWriter writer=null;
		boolean error=false;
		try{
			writer=new ByteStreamWriter(sam, overwrite, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder(BUFFER_BYTES);
			bb.append("@HD\tVN:1.6\tSO:unsorted\n");
			bb.append("@SQ\tSN:chr1\tLN:100\n");
			for(int i=0; i<triplets; i++){
				appendSamRecord(bb, MATCH_PREFIX, i, MATCH_SUFFIX);
				appendSamRecord(bb, SUB_PREFIX, i, SUB_SUFFIX);
				appendSamRecord(bb, DEL_PREFIX, i, DEL_SUFFIX);
				flushIfNeeded(writer, bb);
			}
			flush(writer, bb);
		}finally{
			if(writer!=null){error=writer.poisonAndWait();}
		}
		require(!error, "SAM writer reported an error for "+sam);
		require(new File(sam).length()==expected, "SAM byte count changed for triplets="+triplets+
			"; expected "+expected+" but observed "+new File(sam).length());
		return expected;
	}

	/** Writes one scaled five-column q10/q12 oracle. */
	private long writeQualityOracle(final String fname, final int[][] rows, final long expected){
		ByteStreamWriter writer=null;
		boolean error=false;
		try{
			writer=new ByteStreamWriter(fname, overwrite, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder(512);
			for(final int[] row : rows){
				bb.append(row[0]).tab();
				bb.append(row[1]).tab();
				bb.append(row[2]).tab();
				bb.append(row[3]*1L*triplets).tab();
				bb.append(row[4]*1L*triplets).nl();
			}
			writer.print(bb);
		}finally{
			if(writer!=null){error=writer.poisonAndWait();}
		}
		require(!error, "Oracle writer reported an error for "+fname);
		require(new File(fname).length()==expected, "Oracle byte count changed for "+fname+
			"; expected "+expected+" but observed "+new File(fname).length());
		return expected;
	}

	/** Writes the scaled four-column p-matrix oracle. */
	private long writeCycleOracle(final String fname, final int[][] rows, final long expected){
		ByteStreamWriter writer=null;
		boolean error=false;
		try{
			writer=new ByteStreamWriter(fname, overwrite, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder(256);
			for(final int[] row : rows){
				bb.append(row[0]).tab();
				bb.append(row[1]).tab();
				bb.append(row[2]*1L*triplets).tab();
				bb.append(row[3]*1L*triplets).nl();
			}
			writer.print(bb);
		}finally{
			if(writer!=null){error=writer.poisonAndWait();}
		}
		require(!error, "Oracle writer reported an error for "+fname);
		require(new File(fname).length()==expected, "Oracle byte count changed for "+fname+
			"; expected "+expected+" but observed "+new File(fname).length());
		return expected;
	}

	/** Appends one accepted-literal SAM record with only QNAME made unique. */
	private static void appendSamRecord(final ByteBuilder bb, final String prefix,
			final int ordinal, final String suffix){
		bb.append(prefix).append(ordinal).append(suffix);
	}

	/** Flushes large buffers without changing record boundaries. */
	private static void flushIfNeeded(final ByteStreamWriter writer, final ByteBuilder bb){
		if(bb.length()>=BUFFER_FLUSH_BYTES){flush(writer, bb);}
	}

	/** Flushes buffered bytes to the native writer and reuses the builder. */
	private static void flush(final ByteStreamWriter writer, final ByteBuilder bb){
		if(!bb.isEmpty()){writer.print(bb); bb.clear();}
	}

	/** Returns the exact bytes in the uncompressed scaled SAM. */
	private static long expectedSamBytes(final int count){
		return HEADER_BYTES+
			count*(MATCH_PREFIX.length()+MATCH_SUFFIX.length()+
				SUB_PREFIX.length()+SUB_SUFFIX.length()+
				DEL_PREFIX.length()+DEL_SUFFIX.length())+
			3L*decimalBytesTo(count);
	}

	/** Returns the exact bytes in a scaled five-column q10/q12 TSV. */
	private static long expectedQualityBytes(final int[][] rows, final int count){
		long bytes=0;
		for(final int[] row : rows){
			bytes+=decimalBytes(row[0])+1+
				decimalBytes(row[1])+1+
				decimalBytes(row[2])+1+
				decimalBytes(row[3]*1L*count)+1+
				decimalBytes(row[4]*1L*count)+1;
		}
		return bytes;
	}

	/** Returns the exact bytes in a scaled four-column p-matrix TSV. */
	private static long expectedCycleBytes(final int[][] rows, final int count){
		long bytes=0;
		for(final int[] row : rows){
			bytes+=decimalBytes(row[0])+1+
				decimalBytes(row[1])+1+
				decimalBytes(row[2]*1L*count)+1+
				decimalBytes(row[3]*1L*count)+1;
		}
		return bytes;
	}

	/** Returns decimal digits emitted by appending every value from 0 to count-1. */
	private static long decimalBytesTo(final int count){
		long sum=0, first=0, next=10;
		for(int digits=1; first<count; digits++){
			final long end=Math.min(count, next);
			sum+=(end-first)*digits;
			first=end;
			next*=10;
		}
		return sum;
	}

	/** Returns decimal digits emitted by ByteBuilder.append(long) for a nonnegative value. */
	private static int decimalBytes(final long value){
		require(value>=0, "This helper only writes nonnegative matrix counts: "+value);
		int digits=1;
		for(long x=value; x>=10; x/=10){digits++;}
		return digits;
	}

	/** True if count is one of the finite controller-selected workloads. */
	private static boolean isAllowedTriplets(final int count){
		for(final int allowed : ALLOWED_TRIPLETS){
			if(count==allowed){return true;}
		}
		return false;
	}

	/** Checks the recorded finite byte table, including the N=1 equivalence case. */
	private static boolean checkByteBounds(){
		for(int i=0; i<ALLOWED_TRIPLETS.length; i++){
			final int count=ALLOWED_TRIPLETS[i];
			require(expectedSamBytes(count)==SAM_BYTES[i], "SAM byte bound changed for "+count);
			require(expectedQualityBytes(OLD_Q12_ROWS, count)==OLD_Q12_BYTES[i],
				"old-axis q10/q12 byte bound changed for "+count);
			require(expectedQualityBytes(NEXT_Q12_ROWS, count)==NEXT_Q12_BYTES[i],
				"next-axis q12 byte bound changed for "+count);
			require(expectedCycleBytes(P_ROWS, count)==P_BYTES[i],
				"p-matrix byte bound changed for "+count);
		}
		return true;
	}

	/** True if two paths target the same absolute pathname after normalization. */
	private static boolean sameFile(final File a, final File b){
		return a.getAbsoluteFile().toPath().normalize().equals(b.getAbsoluteFile().toPath().normalize());
	}

	/** Throws one explicit argument or generation error. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Parsed command-line flags. */
	private static final class Args{

		/** Parses strict flag=value arguments. */
		Args(final String[] args){
			final Parser parser=new Parser();
			for(final String arg : args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				require(split.length==2, "Expected flag=value, not "+arg+"\n"+USAGE);
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split[1];
				if(a.equals("out") || a.equals("sam")){sam=b;}
				else if(a.equals("q10")){q10=b;}
				else if(a.equals("q12old")){q12old=b;}
				else if(a.equals("q12next")){q12next=b;}
				else if(a.equals("p")){p=b;}
				else if(a.equals("triplets") || a.equals("count")){triplets=Parse.parseIntKMG(b);}
				else if(a.equals("ow") || a.equals("overwrite")){overwrite=Parse.parseBoolean(b);}
				else if(parser.parse(arg, a, b)){
					// Consumed by BBTools shared parsing.
				}else{
					throw new IllegalArgumentException("Unknown CalcTrueQuality q12 timing-data option: "+arg+"\n"+USAGE);
				}
			}
			require(triplets>0, "Set triplets=1|250k|1m|4m|8m.\n"+USAGE);
		}

		String sam;
		String q10;
		String q12old;
		String q12next;
		String p;
		int triplets=-1;
		boolean overwrite=false;
	}

	/** Usage string with exact uncompressed byte bounds for accepted workloads. */
	private static final String USAGE=
		"Usage: dev/calctruequalityq12timingdata.sh triplets=1|250k|1m|4m|8m "+
		"out=input.sam q10=q10.tsv q12old=q12old.tsv q12next=q12next.tsv p=p.tsv [ow=f]\n"+
		"Exact uncompressed byte bounds by triplets:\n"+
		"1: sam=165 q10=120 q12old=120 q12next=120 p=32\n"+
		"250k: sam=34166712 q10=185 q12old=185 q12next=185 p=65\n"+
		"1m: sam=137666712 q10=198 q12old=198 q12next=198 p=68\n"+
		"4m: sam=560666712 q10=198 q12old=198 q12next=198 p=72\n"+
		"8m: sam=1124666712 q10=209 q12old=209 q12next=209 p=73";

	/** Accepted finite sizes: tiny N=1 plus the four pilot sizes from the timing plan. */
	private static final int[] ALLOWED_TRIPLETS={1, 250000, 1000000, 4000000, 8000000};
	private static final long[] SAM_BYTES={165, 34166712, 137666712, 560666712, 1124666712};
	private static final long[] OLD_Q12_BYTES={120, 185, 198, 198, 209};
	private static final long[] NEXT_Q12_BYTES={120, 185, 198, 198, 209};
	private static final long[] P_BYTES={32, 65, 68, 72, 73};

	/** Fixed accepted header from /tmp/shinobu-jgi/calctruequality-q12-baseline-1/input.sam. */
	private static final int HEADER_BYTES=42;

	/** Accepted non-qname SAM fields. */
	private static final String MATCH_PREFIX="r_match_", SUB_PREFIX="r_sub_", DEL_PREFIX="r_del_";
	private static final String MATCH_SUFFIX="\t0\tchr1\t1\t60\t3=\t*\t0\t0\tACG\t+5?\n";
	private static final String SUB_SUFFIX="\t0\tchr1\t1\t60\t1=1X1=\t*\t0\t0\tACG\t,6@\n";
	private static final String DEL_SUFFIX="\t0\tchr1\t1\t60\t1=1D3=\t*\t0\t0\tACGT\t-7AK\n";

	/** Q10 and buggy q12 rows: pairnum, currentQ, previousQ, sum, bad. */
	private static final int[][] OLD_Q12_ROWS={
		{0, 10, 51, 2, 0},
		{0, 11, 51, 2, 0},
		{0, 12, 51, 2, 1},
		{0, 20, 10, 2, 0},
		{0, 21, 11, 2, 2},
		{0, 22, 12, 2, 1},
		{0, 30, 20, 2, 0},
		{0, 31, 21, 2, 0},
		{0, 32, 22, 2, 0},
		{0, 42, 32, 2, 0},
	};

	/** Fixed q12 rows: pairnum, currentQ, nextQ, sum, bad. */
	private static final int[][] NEXT_Q12_ROWS={
		{0, 10, 20, 2, 0},
		{0, 11, 21, 2, 0},
		{0, 12, 22, 2, 1},
		{0, 20, 30, 2, 0},
		{0, 21, 31, 2, 2},
		{0, 22, 32, 2, 1},
		{0, 30, 51, 2, 0},
		{0, 31, 51, 2, 0},
		{0, 32, 42, 2, 0},
		{0, 42, 51, 2, 0},
	};

	/** P-matrix rows: pairnum, cycle, sum, bad. */
	private static final int[][] P_ROWS={
		{0, 0, 6, 1},
		{0, 1, 6, 3},
		{0, 2, 6, 0},
		{0, 3, 2, 0},
	};

	/** ByteBuilder capacity and flush threshold. */
	private static final int BUFFER_BYTES=1<<20;
	private static final int BUFFER_FLUSH_BYTES=1<<19;

	private final int triplets;
	private final boolean overwrite;
	private final String sam;
	private final String q10;
	private final String q12old;
	private final String q12next;
	private final String p;

	private long samBytes;
	private long q10Bytes;
	private long q12oldBytes;
	private long q12nextBytes;
	private long pBytes;

}
