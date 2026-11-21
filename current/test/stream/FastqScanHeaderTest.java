package test.stream;

import java.io.BufferedOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.security.DigestOutputStream;
import java.security.MessageDigest;
import java.util.Arrays;
import java.util.List;

import fileIO.FileFormat;
import parse.Parse;
import shared.Shared;
import stream.FASTQ;
import stream.FastqScanStreamer;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;

/** Literal raw-header controls for the factory-selected FASTQ scanner.
 * Four fresh JVMs cross single/interleaved input with description trimming off/on.
 * Baseline records the current full-header construction; fixed requires trimming.
 * Empty trimmed prefixes are null, while an untrimmed empty header is an empty String.
 * No processor, writer or numeric-name fallback participates. The default matrix is untimed.
 * Separate write/time modes support a fixed repeated input and actual single-reader timings;
 * their caller owns calibration, repeated samples and input hash verification outside timing.
 * @author Sucrose
 */
public final class FastqScanHeaderTest{

	/** Matrix: phase and new outdir; write: mode=write and new outdir; time: phase, in, trd and passes. */
	public static void main(final String[] args) throws Exception{
		String phase=null, directory=null, input=null, observations=null;
		String mode="matrix";
		int passes=-1;
		boolean child=false;
		Boolean paired=null, trim=null;
		for(final String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			if(pair.length!=2){throw new IllegalArgumentException(arg);}
			if(pair[0].equals("phase")){phase=pair[1];}
			else if(pair[0].equals("mode")){mode=pair[1];}
			else if(pair[0].equals("passes")){passes=Integer.parseInt(pair[1]);}
			else if(pair[0].equals("outdir")){directory=pair[1];}
			else if(pair[0].equals("child")){child=Parse.parseBoolean(pair[1]);}
			else if(pair[0].equals("in")){input=pair[1];}
			else if(pair[0].equals("observations")){observations=pair[1];}
			else if(pair[0].equals("paired")){paired=Parse.parseBoolean(pair[1]);}
			else if(pair[0].equals("trd")){trim=Parse.parseBoolean(pair[1]);}
			else{throw new IllegalArgumentException(arg);}
		}
		if(!mode.equals("matrix") && !mode.equals("write") && !mode.equals("time")){
			throw new IllegalArgumentException("Set mode=matrix|write|time.");
		}
		if(!mode.equals("write") && !"baseline".equals(phase) && !"fixed".equals(phase)){
			throw new IllegalArgumentException("Set phase=baseline|fixed.");
		}
		boolean assertions=false;
		assert(assertions=true);
		check(assertions, "Reader controls require assertions enabled in every JVM.");
		check(Charset.defaultCharset().equals(StandardCharsets.UTF_8),
			"The literal UTF-8 input tests existing default decoding; require -Dfile.encoding=UTF-8.");
		if(mode.equals("write")){
			if(child || phase!=null || directory==null || input!=null || observations!=null || paired!=null || trim!=null || passes!=-1){
				throw new IllegalArgumentException("Write mode requires only mode=write and a new outdir.");
			}
			writeTimingInput(Files.createDirectory(Paths.get(directory)).toAbsolutePath());
		}else if(mode.equals("time")){
			if(child || directory!=null || input==null || observations!=null || paired!=null || trim==null || passes<1){
				throw new IllegalArgumentException("Time mode requires phase, in, trd and positive passes; single reads only.");
			}
			timeReader(Paths.get(input), phase, trim, passes);
		}else if(child){
			if(directory!=null || input==null || observations==null || paired==null || trim==null || passes!=-1){
				throw new IllegalArgumentException("Child requires in, observations, paired and trd, without outdir.");
			}
			readCase(Paths.get(input), Paths.get(observations), phase, paired, trim);
		}else{
			if(directory==null || input!=null || observations!=null || paired!=null || trim!=null || passes!=-1){
				throw new IllegalArgumentException("Parent requires only phase and a new outdir; run from the BBTools root.");
			}
			runMatrix(Files.createDirectory(Paths.get(directory)).toAbsolutePath(), phase);
		}
	}

	/** Saves one literal input and each exact oracle before starting its fresh child JVM. */
	private static void runMatrix(final Path root, final String phase) throws Exception{
		final Path input=root.resolve("headers.fq");
		Files.write(input, inputBlock());
		final Path launcher=Paths.get("dev/fastqscanheadertest.sh").toAbsolutePath();
		if(!Files.isRegularFile(launcher)){throw new IllegalArgumentException("Run from the BBTools root: "+launcher);}
		for(final boolean paired : new boolean[]{false, true}){
			for(final boolean trim : new boolean[]{false, true}){
				final String label=(paired ? "interleaved" : "single")+"-trd"+trim;
				final Path actual=root.resolve(label+".observed.tsv");
				final byte[] expected=expectedRows(phase, paired, trim).getBytes(StandardCharsets.UTF_8);
				Files.write(root.resolve(label+".expected.tsv"), expected);
				final List<String> command=Arrays.asList("bash", launcher.toString(), "child=t", "phase="+phase,
					"in="+input, "observations="+actual, "paired="+paired, "trd="+trim);
				Files.write(root.resolve(label+".args"), command, StandardCharsets.UTF_8);
				final Path stdout=root.resolve(label+".stdout"), stderr=root.resolve(label+".stderr");
				final Process process=new ProcessBuilder(command).redirectOutput(stdout.toFile()).redirectError(stderr.toFile()).start();
				final int code=process.waitFor();
				Files.write(root.resolve(label+".exit"), (code+"\n").getBytes(StandardCharsets.US_ASCII));
				check(code==0, label+": reader child failed; inspect saved observations and stderr.");
				final byte[] marker=(caseMarker(phase, paired, trim)+"\n").getBytes(StandardCharsets.US_ASCII);
				check(Arrays.equals(marker, Files.readAllBytes(stdout)), label+": missing exact child completion/count marker.");
				check(Arrays.equals(expected, Files.readAllBytes(actual)), label+": raw reader observations differ from literal expectations.");
			}
		}
		check(checks==PARENT_CHECKS-1, "Parent check count changed: "+checks+"; review the frozen matrix before accepting it.");
		System.out.println("FASTQ_SCAN_HEADER_PASS phase="+phase+" checks="+(checks+4*CHILD_CHECKS)+
			" parentChecks="+checks+" childChecks="+(4*CHILD_CHECKS)+" children=4");
	}

	/** Configures the actual factory path, drains raw Reads and retains evidence before each assertion. */
	private static void readCase(final Path input, final Path observations, final String phase,
		final boolean paired, final boolean trim) throws Exception{
		final Streamer reader=makeReader(input, paired, trim);
		final String[] names=phase.equals("fixed") && trim ? TRIMMED : HEADERS;
		int entries=0, records=0, bases=0;
		try(PrintWriter out=new PrintWriter(Files.newBufferedWriter(observations, StandardCharsets.UTF_8,
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE))){
			out.print(COLUMNS);
			try{
				reader.start();
				for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
					for(final Read read : batch.list){
						verifyRead(out, read, records, paired, names, paired ? read.mate : null);
						entries++; records++; bases+=read.length();
						if(paired){
							verifyRead(out, read.mate, records, true, names, read);
							records++; bases+=read.mate.length();
						}
					}
				}
				check(entries==(paired ? 4 : 8), "Exactly eight input records must become eight singles or four paired roots; roots="+entries);
				check(records==8, "Every literal header must be observed exactly once; records="+records);
				check(bases==32, "The eight four-base payloads must retain all32 bases; bases="+bases);
				check(!reader.hasMore(), "The fixture must consume the reader's terminal marker.");
				check(!reader.errorState(), "The complete literal FASTQ must not leave a reader error.");
				check(((FastqScanStreamer)reader).corruption()==null, "The complete LF-terminated fixture must have no structural corruption.");
			}finally{reader.close();}
			if(out.checkError()){throw new AssertionError("Could not retain raw reader observations.");}
		}
		check(checks==CHILD_CHECKS-1, "Child check count changed: "+checks+"; review the frozen record contract.");
		System.out.println(caseMarker(phase, paired, trim));
	}

	/** Forces the same actual scanner route for finite controls and each fresh timing pass. */
	private static Streamer makeReader(final Path input, final boolean paired, final boolean trim){
		Shared.setThreads(1);
		Shared.SIMD=false;
		Shared.TRIM_READ_DESCRIPTION=trim;
		FASTQ.FORCE_INTERLEAVED=paired;
		FASTQ.TEST_INTERLEAVED=false;
		FASTQ.ASCII_OFFSET=33;
		FASTQ.DETECT_QUALITY=false;
		FASTQ.SET_QIN=true;
		Read.FIX_HEADER=false;
		final FileFormat format=FileFormat.testInput(input.toString(), FileFormat.FASTQ, null, false, true);
		check(format.fastq() && format.interleaved()==paired, "The input descriptor must retain the requested FASTQ/interleaving mode.");
		final Streamer reader=StreamerFactory.makeStreamer(format, 0, true, -1, false, true, 1);
		check(reader instanceof FastqScanStreamer, "This control must reach FastqScanStreamer, not another decoder.");
		check(reader.paired()==paired, "Reader construction must retain the requested interleaving mode.");
		return reader;
	}

	/** Writes exactly65,536 literal eight-read blocks and hashes bytes during their one write. */
	private static void writeTimingInput(final Path root) throws Exception{
		final byte[] block=inputBlock();
		final MessageDigest digest=MessageDigest.getInstance("SHA-256");
		final Path input=root.resolve("headers.fq");
		try(OutputStream out=new BufferedOutputStream(new DigestOutputStream(Files.newOutputStream(input,
			StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE), digest))){
			for(int i=0; i<TIMING_BLOCKS; i++){out.write(block);}
		}
		final long bytes=(long)block.length*TIMING_BLOCKS;
		if(Files.size(input)!=bytes){throw new AssertionError("Fixed timing input has an unexpected byte count.");}
		final String spec="encoding\tUTF-8\nblocks\t"+TIMING_BLOCKS+"\nreads\t"+TIMING_READS+
			"\nbases\t"+TIMING_BASES+"\nbytes\t"+bytes+"\nsha80\t"+sha80(digest.digest())+"\n";
		Files.write(root.resolve("input-spec.tsv"), spec.getBytes(StandardCharsets.US_ASCII));
		System.out.println("FASTQ_SCAN_HEADER_INPUT_PASS blocks="+TIMING_BLOCKS+" reads="+TIMING_READS+
			" bases="+TIMING_BASES+" bytes="+bytes);
	}

	/** Times repeated real file readers, including construction, queues, Read allocation and a result digest.
	 * One literal expected-pass digest is computed before the pass loop; no input hash is computed here.
	 * Native whole-JVM timing also includes this fixed oracle work and JVM startup. The caller calibrates
	 * passes until actual compute exceeds ten seconds, then keeps the same passes before/after repair.
	 */
	private static void timeReader(final Path input, final String phase, final boolean trim, final int passes){
		final long expectedPass=expectedPassDigest(phase, trim);
		long digest=HASH_SEED, expectedDigest=HASH_SEED, totalReads=0, totalBases=0;
		for(int pass=0; pass<passes; pass++){
			final Streamer reader=makeReader(input, false, trim);
			long passDigest=HASH_SEED, reads=0, bases=0;
			try{
				reader.start();
				for(ListNum<Read> batch=reader.nextList(); batch!=null; batch=reader.nextList()){
					for(final Read read : batch.list){
						passDigest=hashRead(passDigest, read.id, read.bases, read.quality);
						reads++; bases+=read.length();
					}
				}
				if(reader.hasMore() || reader.errorState() || ((FastqScanStreamer)reader).corruption()!=null){
					throw new AssertionError("Timing pass "+pass+" did not drain a complete valid input.");
				}
			}finally{reader.close();}
			if(reads!=TIMING_READS || bases!=TIMING_BASES || passDigest!=expectedPass){
				throw new AssertionError("Timing pass "+pass+" differs: reads="+reads+"/"+TIMING_READS+
					", bases="+bases+"/"+TIMING_BASES+", digest="+Long.toHexString(passDigest)+"/"+Long.toHexString(expectedPass));
			}
			totalReads+=reads; totalBases+=bases;
			digest=(digest^passDigest)*HASH_MULTIPLIER;
			expectedDigest=(expectedDigest^expectedPass)*HASH_MULTIPLIER;
		}
		System.out.println("FASTQ_SCAN_HEADER_TIME_PASS phase="+phase+" trd="+trim+" passes="+passes+
			" readsPerPass="+TIMING_READS+" basesPerPass="+TIMING_BASES+" reads="+totalReads+" bases="+totalBases+
			" expectedDigest="+Long.toHexString(expectedDigest)+" digest="+Long.toHexString(digest)+
			" passDigest="+Long.toHexString(expectedPass));
	}

	/** Calculates the expected ordered digest from literals, without opening an input or calling a reader. */
	private static long expectedPassDigest(final String phase, final boolean trim){
		final String[] names=phase.equals("fixed") && trim ? TRIMMED : HEADERS;
		final byte[][] sequences=new byte[8][];
		for(int i=0; i<8; i++){sequences[i]=SEQUENCES[i].getBytes(StandardCharsets.US_ASCII);}
		long hash=HASH_SEED;
		for(int block=0; block<TIMING_BLOCKS; block++){
			for(int i=0; i<8; i++){hash=hashRead(hash, names[i], sequences[i], QUALITIES);}
		}
		return hash;
	}

	/** FNV-style64-bit result checksum; length tags distinguish null/empty IDs and field boundaries.
	 * Hashes decoded UTF-16 ID units and literal base/quality bytes; not a cryptographic input hash.
	 */
	private static long hashRead(long hash, final String id, final byte[] bases, final byte[] qualities){
		hash=(hash^(id==null ? -1 : id.length()))*HASH_MULTIPLIER;
		if(id!=null){for(int i=0; i<id.length(); i++){hash=(hash^id.charAt(i))*HASH_MULTIPLIER;}}
		hash=(hash^(bases==null ? -1 : bases.length))*HASH_MULTIPLIER;
		if(bases!=null){for(final byte base : bases){hash=(hash^(base&255))*HASH_MULTIPLIER;}}
		hash=(hash^(qualities==null ? -1 : qualities.length))*HASH_MULTIPLIER;
		if(qualities!=null){for(final byte quality : qualities){hash=(hash^(quality&255))*HASH_MULTIPLIER;}}
		return hash;
	}

	/** Formats only the final80 bits of the input SHA-256; the full hash is never printed. */
	private static String sha80(final byte[] digest){
		final StringBuilder text=new StringBuilder(20);
		for(int i=digest.length-10; i<digest.length; i++){
			text.append(Character.forDigit((digest[i]>>>4)&15, 16)).append(Character.forDigit(digest[i]&15, 16));
		}
		return text.toString();
	}

	/** The same eight literal UTF-8 records feed the finite matrix and repeated timing input. */
	private static byte[] inputBlock(){
		final StringBuilder text=new StringBuilder();
		for(int i=0; i<HEADERS.length; i++){
			text.append('@').append(HEADERS[i]).append('\n').append(SEQUENCES[i]).append("\n+\nIIII\n");
		}
		return text.toString().getBytes(StandardCharsets.UTF_8);
	}

	/** Checks six independent raw-read properties; never derives expected names from reader code. */
	private static void verifyRead(final PrintWriter out, final Read read, final int index,
		final boolean paired, final String[] names, final Read mate){
		if(read==null || index>=8){throw new AssertionError("Missing mate or unexpected extra read at input position "+index);}
		out.print(row(index, read.id, new String(read.bases, StandardCharsets.US_ASCII), Arrays.toString(read.quality),
			read.numericID, read.pairnum(), read.mate==null ? "<null>" : read.mate.numericID+":"+read.mate.pairnum(),
			read.mate!=null && read.mate.mate==read));
		out.flush();
		if(out.checkError()){throw new AssertionError("Could not save raw record "+index);}
		final String label="Input record "+index+": ";
		check(names[index]==null ? read.id==null : names[index].equals(read.id),
			label+"raw ID differs; expected="+escape(names[index])+", actual="+escape(read.id));
		check(Arrays.equals(SEQUENCES[index].getBytes(StandardCharsets.US_ASCII), read.bases), label+"sequence changed.");
		check(Arrays.equals(QUALITIES, read.quality), label+"numeric qualities must remain [40,40,40,40].");
		check(read.numericID==(paired ? index/2 : index), label+"input-position numeric ID changed: "+read.numericID);
		check(read.pairnum()==(paired ? index%2 : 0), label+"pair number changed: "+read.pairnum());
		check(paired ? mate!=null && read.mate==mate && mate.mate==read : read.mate==null,
			label+"single reads must have no mate; paired reads must retain reciprocal links.");
	}

	/** Literal expected records for the two explicitly separate behavior phases. */
	private static String expectedRows(final String phase, final boolean paired, final boolean trim){
		final String[] names=phase.equals("fixed") && trim ? TRIMMED : HEADERS;
		final StringBuilder text=new StringBuilder(COLUMNS);
		for(int i=0; i<8; i++){
			text.append(row(i, names[i], SEQUENCES[i], "[40, 40, 40, 40]", paired ? i/2 : i,
				paired ? i%2 : 0, paired ? (i/2)+":"+(1-i%2) : "<null>", paired));
		}
		return text.toString();
	}

	/** Serializes evidence only; quoted empty strings and null remain distinguishable. */
	private static String row(final int index, final String id, final String bases, final String qualities,
		final long numericID, final int pairnum, final String mate, final boolean reciprocal){
		return index+"\t"+escape(id)+"\t"+bases+"\t"+qualities+"\t"+numericID+"\t"+pairnum+"\t"+mate+"\t"+reciprocal+"\n";
	}

	/** Escapes evidence delimiters without changing the ID supplied to the reader oracle. */
	private static String escape(final String text){
		return text==null ? "<null>" : "\""+text.replace("\\", "\\\\").replace("\"", "\\\"")
			.replace("\t", "\\t").replace("\n", "\\n")+"\"";
	}

	/** Exact marker checked by the parent after the child has drained and closed its input. */
	private static String caseMarker(final String phase, final boolean paired, final boolean trim){
		return "FASTQ_SCAN_HEADER_CASE_PASS phase="+phase+" paired="+paired+" trd="+trim+" checks="+CHILD_CHECKS;
	}

	/** Literal oracles remain active even if a later launcher drops assertions; main rejects that setup. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	/** Headers after exactly one physical FASTQ marker; the fourth input therefore begins with @@. */
	private static final String[] HEADERS={"plain", "space desc", "tab\tdesc", "@embedded desc", "\u00e9 desc", "", " leading", "\tleading"};
	/** Intended trimmed prefixes are literals, including three null raw IDs. */
	private static final String[] TRIMMED={"plain", "space", "tab", "@embedded", "\u00e9", null, null, null};
	/** Distinct literal sequences expose dropped, duplicated or reordered records. */
	private static final String[] SEQUENCES={"ACGT", "TGCA", "AAAA", "CCCC", "GGGG", "TTTT", "AGCT", "TCGA"};
	private static final byte[] QUALITIES={40, 40, 40, 40};
	private static final String COLUMNS="record\tid\tbases\tqualities\tnumericID\tpairnum\tmate\treciprocal\n";
	//Child:2 JVM preconditions+3 route/mode+8*6 read properties+6 terminal/totals+1 sentry=60.
	//Parent:2 JVM preconditions+4*3 child artifact checks+1 sentry=15; total255.
	private static final int CHILD_CHECKS=60, PARENT_CHECKS=15;
	/** Fixed workload size; pass count alone is calibrated by the runner. */
	private static final int TIMING_BLOCKS=65536;
	private static final long TIMING_READS=8L*TIMING_BLOCKS, TIMING_BASES=32L*TIMING_BLOCKS;
	private static final long HASH_SEED=0xcbf29ce484222325L, HASH_MULTIPLIER=0x100000001b3L;
	private static int checks;

	private FastqScanHeaderTest(){}
}
