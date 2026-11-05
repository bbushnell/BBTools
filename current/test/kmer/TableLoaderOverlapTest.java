package test.kmer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;

import fileIO.FileFormat;
import jgi.BBMerge;
import kmer.HashArray1D;
import kmer.TableLoaderLockFree;
import stream.FASTQ;
import stream.FastqScanStreamer;
import stream.FastqStreamer;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Separates paired-reference correction from optional name retention.
 * One 35-base pair has a single low-quality error opposite a correct Q40 base.
 * Real-reader and strict-overlap guards qualify the input before loader assertions.
 * Ten frozen windows and independent encoding define the expected stored keys.
 * Run one case per fresh JVM under either assertion mode; this is not a timing test.
 * @author Sucrose
 */
public final class TableLoaderOverlapTest{

	/** Accepts a case and a directory in which the tiny FASTQ is retained without overwrite. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path dir=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("dir=")){dir=Paths.get(arg.substring(4));}
			else{throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		fixture(mode!=null && dir!=null, mode, "Expected case= and dir= arguments");
		final boolean ecc, keepNames;
		switch(mode){
			case "ecc-no-names": ecc=true; keepNames=false; break;
			case "ecc-names-control": ecc=true; keepNames=true; break;
			case "no-ecc-no-names-control": ecc=false; keepNames=false; break;
			case "no-ecc-names-control": ecc=false; keepNames=true; break;
			default: throw new IllegalArgumentException("Unknown overlap case: "+mode);
		}
		final long[] keys=expectedKeys(mode);
		Files.createDirectories(dir);
		final Path input=dir.resolve("paired-reference.fq");
		final byte[] quality1=new byte[35], quality2=new byte[35];
		Arrays.fill(quality1, (byte)'I');
		Arrays.fill(quality2, (byte)'I');
		quality1[17]='#';//Phred33 Q2; all other calls are Q40.
		final String fastq="@alpha\n"+ERROR+"\n+\n"+new String(quality1, StandardCharsets.US_ASCII)+
			"\n@unrelated\n"+REVERSE+"\n+\n"+new String(quality2, StandardCharsets.US_ASCII)+"\n";
		Files.write(input, fastq.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);

		final boolean oldForce=FASTQ.FORCE_INTERLEAVED, oldTest=FASTQ.TEST_INTERLEAVED;
		final boolean oldDetect=FASTQ.DETECT_QUALITY, oldProgress=TableLoaderLockFree.DISPLAY_PROGRESS;
		final boolean oldAmbiguous=TableLoaderLockFree.REPLICATE_AMBIGUOUS;
		final byte oldOffset=FASTQ.ASCII_OFFSET;
		final int oldReaderThreads=FastqStreamer.DEFAULT_THREADS;
		try{
			FASTQ.FORCE_INTERLEAVED=true;
			FASTQ.TEST_INTERLEAVED=false;
			FASTQ.ASCII_OFFSET=33;
			FASTQ.DETECT_QUALITY=false;
			FastqStreamer.DEFAULT_THREADS=1;//Factory selects the scan reader regardless of global worker count.
			TableLoaderLockFree.DISPLAY_PROGRESS=false;
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=false;
			qualify(input, mode);
			load(input, keys, ecc, keepNames, mode);
		}finally{
			FASTQ.FORCE_INTERLEAVED=oldForce;
			FASTQ.TEST_INTERLEAVED=oldTest;
			FASTQ.ASCII_OFFSET=oldOffset;
			FASTQ.DETECT_QUALITY=oldDetect;
			FastqStreamer.DEFAULT_THREADS=oldReaderThreads;
			TableLoaderLockFree.DISPLAY_PROGRESS=oldProgress;
			TableLoaderLockFree.REPLICATE_AMBIGUOUS=oldAmbiguous;
		}
		System.out.println("PASS loader-overlap case="+mode+" keys="+(ecc ? 5 : 10));
	}

	/** Requires successful strict overlap with and without correction on separately read pairs. */
	private static void qualify(final Path input, final String mode){
		final Read unchanged=readPair(input, mode);
		fixture(BBMerge.findOverlapStrict(unchanged, unchanged.mate, false)==35, mode,
			"Strict overlap without correction must find insert35 before testing the loader");
		fixture(bases(unchanged, ERROR) && bases(unchanged.mate, REVERSE), mode,
			"Overlap without correction must preserve both input sequences and mate orientation");
		final Read corrected=readPair(input, mode);
		fixture(BBMerge.findOverlapStrict(corrected, corrected.mate, true)==35, mode,
			"Strict correction must find insert35 before testing the loader");
		fixture(bases(corrected, CORRECT) && bases(corrected.mate, REVERSE), mode,
			"Strict correction must repair exactly the R1 error and restore R2 orientation");
	}

	/** Uses the same descriptor/factory call as the loader; pairing never depends on header agreement. */
	private static Read readPair(final Path input, final String mode){
		final FileFormat ff=FileFormat.testInput(input.toString(), FileFormat.FASTA, null, false, true);
		fixture(ff.fastq() && ff.interleaved(), mode, "The reference descriptor must select interleaved FASTQ");
		final Streamer reader=StreamerFactory.getReadInputStream(-1, false, ff, null, -1);
		fixture(reader instanceof FastqScanStreamer && reader.paired(), mode,
			"Factory must select the paired scan reader, independently of headers");
		Read first=null;
		int roots=0;
		try{
			reader.start();
			for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
				fixture(batch.list!=null, mode, "A nonterminal batch must have a list");
				for(Read r : batch.list){
					roots++;
					fixture(roots==1 && r!=null, mode, "Input must contain exactly one nonnull pair root");
					first=r;
				}
			}
		}finally{reader.close();}
		fixture(!reader.errorState() && roots==1 && first!=null && first.mate!=null, mode,
			"Successful real input must produce exactly one pair");
		fixture(first.mate.mate==first && first.pairnum()==0 && first.mate.pairnum()==1, mode,
			"Input must link reciprocal first/second mates");
		fixture(bases(first, ERROR) && bases(first.mate, REVERSE), mode, "Reader must preserve the two literal sequences");
		fixture(first.quality!=null && first.mate.quality!=null && first.quality.length==35
			&& first.mate.quality.length==35, mode, "Each 35-base mate needs 35 decoded qualities");
		for(int i=0; i<35; i++){
			fixture(first.quality[i]==(i==17 ? 2 : 40) && first.mate.quality[i]==40, mode,
				"Real Phred33 decoding must preserve the intended quality difference at index "+i);
		}
		return first;
	}

	/** Checks the loader against frozen windows; the erroneous-key witness names the intended failure. */
	private static void load(final Path input, final long[] keys, final boolean ecc,
			final boolean keepNames, final String mode) throws ReflectiveOperationException{
		final HashArray1D[] tables=new HashArray1D[7];
		for(int i=0; i<tables.length; i++){tables[i]=new HashArray1D(new int[] {17}, -1L);}
		fixture(storedKeys(tables, mode)==0, mode, "Loader requires seven empty tiny shards");
		final TableLoaderLockFree loader=new TableLoaderLockFree(tables, 31, 0, 0, 0, 0, true, false);
		loader.setRefSkip(0);
		loader.storeMode(TableLoaderLockFree.SET_IF_NOT_PRESENT);
		final long added=loader.processData(new String[] {input.toString()}, null, keepNames, false, ecc);
		final long stored=storedKeys(tables, mode);
		fixture(!loader.errorState && statistic(loader, "refReads")==2 && statistic(loader, "refBases")==70
			&& statistic(loader, "refKmers")==10, mode, "Loader must count two reads,70 bases and10 full windows without errors");
		for(int i=0; i<5; i++){
			fixture(value(tables, keys[i])==1, mode, "Correct mate must supply every correct key with ID1: "+keys[i]);
		}
		if(ecc && value(tables, keys[5])!=-1){
			throw new AssertionError("KMR123 case="+mode+": Requested correction must remove the erroneous literal key");
		}
		for(int i=5; i<10; i++){
			fixture(value(tables, keys[i])==(ecc ? -1 : 1), mode,
				"Every erroneous window must be absent with ECC and retain ID1 without it: "+keys[i]);
		}
		final int expected=ecc ? 5 : 10;
		fixture(stored==expected && added==expected, mode,
			"Exact stored set and creation return must agree: expected="+expected+", stored="+stored+", added="+added);
	}

	/** Independently encodes the ten fixed strings and rejects any duplicate canonical key. */
	private static long[] expectedKeys(final String mode){
		fixture(CORRECT.length()==35 && ERROR.length()==35 && REVERSE.length()==35, mode,
			"Strict overlap requires35 bases in each literal");
		final long[] keys=new long[10];
		for(int i=0; i<10; i++){
			final String window=WINDOWS[i];
			fixture(window.length()==31 && window.equals((i<5 ? CORRECT : ERROR).substring(i%5, i%5+31)), mode,
				"Frozen expected window must match its stated literal position: "+i);
			long forward=0, reverse=0;
			for(int j=0; j<31; j++){
				forward=(forward<<2)|base(window.charAt(j));
				reverse=(reverse<<2)|(3-base(window.charAt(30-j)));
			}
			keys[i]=Math.max(forward, reverse)|(1L<<62);
			for(int j=0; j<i; j++){fixture(keys[i]!=keys[j], mode, "All ten literal canonical windows must be distinct");}
		}
		return keys;
	}

	/** Defines the test's two-bit alphabet independently of production tables/helpers. */
	private static int base(final char c){
		switch(c){
			case 'A': return 0;
			case 'C': return 1;
			case 'G': return 2;
			case 'T': return 3;
			default: throw new AssertionError("Fixture literal contains non-ACGT base: "+c);
		}
	}

	/** Counts actual occupied cells after checking fixed capacity, shard placement and no victims. */
	private static long storedKeys(final HashArray1D[] tables, final String mode){
		long total=0;
		for(int shard=0; shard<7; shard++){
			final HashArray1D table=tables[shard];
			fixture(table.arrayLength()==77 && table.victims().size()==0, mode, "Tiny shards must not grow or overflow");
			long occupied=0;
			for(long key : table.array()){
				if(key>=0){
					fixture(key%7==shard, mode, "Stored key must be in its modulo-seven shard");
					occupied++;
				}
			}
			fixture(occupied==table.size(), mode, "Actual occupancy must match reported primary size");
			total+=occupied;
		}
		return total;
	}

	/** Looks up the fixed seven-way shard without invoking loader key construction. */
	private static int value(final HashArray1D[] tables, final long key){return tables[(int)(key%7)].getValue(key);}

	/** Compares literal bases, including the reader/merger's restored mate orientation. */
	private static boolean bases(final Read r, final String literal){return Arrays.equals(r.bases, literal.getBytes(StandardCharsets.US_ASCII));}

	/** Reads existing aggregate counters without adding production visibility or hooks. */
	private static long statistic(final TableLoaderLockFree loader, final String name) throws ReflectiveOperationException{
		final Field field=TableLoaderLockFree.class.getDeclaredField(name);
		field.setAccessible(true);
		return field.getLong(loader);
	}

	/** Qualification failures remain distinct from the expected KMR123 witness under -ea and -da. */
	private static void fixture(final boolean ok, final String mode, final String message){
		if(!ok){throw new AssertionError("Fixture qualification case="+mode+": "+message);}
	}

	private static final String CORRECT="ACGTTGCAAGTCGATCGGATACCTGAGCTTAGGCA";
	private static final String ERROR="ACGTTGCAAGTCGATCGAATACCTGAGCTTAGGCA";
	private static final String REVERSE="TGCCTAAGCTCAGGTATCCGATCGACTTGCAACGT";
	/** First five belong to the corrected molecule; final five to the erroneous R1. */
	private static final String[] WINDOWS={
		"ACGTTGCAAGTCGATCGGATACCTGAGCTTA",
		"CGTTGCAAGTCGATCGGATACCTGAGCTTAG",
		"GTTGCAAGTCGATCGGATACCTGAGCTTAGG",
		"TTGCAAGTCGATCGGATACCTGAGCTTAGGC",
		"TGCAAGTCGATCGGATACCTGAGCTTAGGCA",
		"ACGTTGCAAGTCGATCGAATACCTGAGCTTA",
		"CGTTGCAAGTCGATCGAATACCTGAGCTTAG",
		"GTTGCAAGTCGATCGAATACCTGAGCTTAGG",
		"TTGCAAGTCGATCGAATACCTGAGCTTAGGC",
		"TGCAAGTCGATCGAATACCTGAGCTTAGGCA"
	};
}
