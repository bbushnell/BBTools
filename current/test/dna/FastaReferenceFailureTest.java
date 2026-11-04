package test.dna;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipException;

import dna.ChromosomeArray;
import dna.FastaToChromArrays2;
import fileIO.ByteFile1;
import fileIO.ReadWrite;
import shared.Shared;

/** Checks the live reference converter with literal input, corrupt gzip and failed metadata output.
 * Valid/CRC cases retain chromosome and metadata output in memory. The metadata
 * case writes only tiny metadata files under its fresh work directory, with info.txt
 * linked to /dev/full. No case writes chromosome objects or reference index data.
 * Restore cases use four-base input or catchable failures before input/writer
 * construction, checking temporary write mode separately from persistent options.
 * @author Sucrose
 */
public final class FastaReferenceFailureTest{

	/** Runs case=valid|crc|metadata|restore-control|restore-missing|restore-parse in a fresh JVM/work directory. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		final boolean restore="restore-control".equals(mode) || "restore-missing".equals(mode) || "restore-parse".equals(mode);
		require(("valid".equals(mode) || "crc".equals(mode) || "metadata".equals(mode) || restore) && work!=null,
				"Require a supported case= and work=<fresh-directory>");
		Files.createDirectory(work);
		Shared.setThreads(1);
		ReadWrite.RAWMODE=false;
		ReadWrite.USE_UNPIGZ=ReadWrite.USE_GUNZIP=ReadWrite.USE_UNBGZIP=false;
		ReadWrite.ALLOW_NATIVE_BGZF=false;
		if(restore){
			writeMode(work, mode);
			System.out.println("PASS "+mode);
			return;
		}
		if(mode.equals("metadata")){
			metadata(work);
			System.out.println("PASS metadata");
			return;
		}
		final byte[] raw=payload();
		require(raw.length>ByteFile1.bufferlen, "Reference must exceed the ByteFile1 buffer");
		final byte[] compressed=gzip(raw);
		verifyGzip(compressed, raw, false);
		if(mode.equals("crc")){
			compressed[compressed.length-8]^=1;
			verifyGzip(compressed, raw, true);
		}
		final Path input=work.resolve("reference.fa.gz");
		Files.write(input, compressed);
		final String[] convert={"in="+input, "build=17", "path="+work,
				"nodisk=t", "writechroms=f", "retain=t", "genscaffoldinfo=t",
				"startpad=0", "midpad=0", "stoppad=0", "minscaf=1", "mergescaffolds=t"};
		ArrayList<ChromosomeArray> result=null;
		RuntimeException failure=null;
		try{result=FastaToChromArrays2.main2(convert);}
		catch(RuntimeException problem){failure=problem;}
		System.out.println("case="+mode+" failure="+failure+" returnedChromosomes="+
				(result==null ? -1 : result.size()));
		if(mode.equals("crc")){
			require(failure!=null, "Corrupt reference returned normally instead of rejecting partial input");
		}else{
			if(failure!=null){throw failure;}
			require(result!=null && result.size()==1, "Expected one retained chromosome");
			final ChromosomeArray ca=result.get(0);
			require(ca.minIndex==0 && ca.maxIndex==BASES-1 && ca.array.length==BASES,
					"Unexpected retained reference bounds or capacity");
			for(int i=0; i<BASES; i++){
				if(ca.get(i)!=MOTIF[i%MOTIF.length]){
					throw new AssertionError("Reference base differs at "+i);
				}
			}
			require(FastaToChromArrays2.INFO_LIST.contains("1\t1\t1\t100000\t100000\t0\t0\t0"),
					"Literal chromosome metadata differs: "+FastaToChromArrays2.INFO_LIST);
			require(FastaToChromArrays2.SCAF_LIST.contains("1\t1\t0\t100000\tone"),
					"Literal scaffold metadata differs: "+FastaToChromArrays2.SCAF_LIST);
			System.out.println("validatedBases="+BASES+" bounds=[0,"+BASES+")");
		}
		require(!Files.exists(work.resolve("ref")), "Memory-only conversion wrote a reference tree");
		System.out.println("PASS "+mode);
	}

	/** Exercises both entry states with default mode and an explicit option which flips that state. */
	private static void writeMode(Path work, String mode) throws Exception{
		final Path input=work.resolve("tiny.fa"), missing=work.resolve("absent.fa");
		Files.write(input, ">one\nACGT\n".getBytes(StandardCharsets.US_ASCII));
		require(!Files.exists(missing), "Missing-input fault path must not exist");
		boolean pass=true;
		for(boolean initial : new boolean[]{false, true}){
			pass&=writeModeCall(work, mode, input, missing, initial, false);
			pass&=writeModeCall(work, mode, input, missing, initial, true);
		}
		require(pass, "Write mode must restore its entry state while unrelated parsed options persist");
	}

	/** Failures occur in build parsing or input-path validation, before conversion can start workers. */
	private static boolean writeModeCall(Path work, String mode, Path input, Path missing, boolean initial, boolean explicit){
		FastaToChromArrays2.WRITE_IN_THREAD=initial;
		FastaToChromArrays2.MIN_SCAFFOLD=1;
		FastaToChromArrays2.MERGE_SCAFFOLDS=true;
		final ArrayList<String> args=new ArrayList<String>();
		args.add("path="+work);
		args.add("in="+(mode.equals("restore-missing") ? missing : input));
		args.add("nodisk=t");
		args.add("writechroms=f");
		args.add("retain=t");
		args.add("genscaffoldinfo=t");
		args.add("startpad=0");
		args.add("midpad=0");
		args.add("stoppad=0");
		args.add("minscaf=2");
		args.add("mergescaffolds=f");
		if(explicit){args.add("writeinthread="+!initial);}
		args.add("build="+(mode.equals("restore-parse") ? "bad-build" : "17"));
		ArrayList<ChromosomeArray> result=null;
		RuntimeException failure=null;
		try{result=FastaToChromArrays2.main2(args.toArray(new String[args.size()]));}
		catch(RuntimeException problem){failure=problem;}
		final boolean success=mode.equals("restore-control");
		final boolean expectedOutcome=success ? failure==null && tinyReference(result) :
				mode.equals("restore-parse") ? failure instanceof NumberFormatException &&
				failure.getMessage()!=null && failure.getMessage().contains("bad-build") :
				failure!=null && failure.getMessage()!=null && failure.getMessage().startsWith("Not a valid file:");
		final boolean restored=FastaToChromArrays2.WRITE_IN_THREAD==initial;
		final boolean persistent=FastaToChromArrays2.MIN_SCAFFOLD==2 && !FastaToChromArrays2.MERGE_SCAFFOLDS;
		final boolean noReferenceTree=!Files.exists(work.resolve("ref"));
		final boolean pass=expectedOutcome && restored && persistent && noReferenceTree;
		System.out.println("case="+mode+" initial="+initial+" explicit="+explicit+
				" requested="+(explicit ? Boolean.toString(!initial) : "default")+
				" exception="+(failure==null ? "none" : failure.getClass().getSimpleName())+
				" outcomeCorrect="+expectedOutcome+" writeMode="+FastaToChromArrays2.WRITE_IN_THREAD+
				" restored="+restored+" minscaf="+FastaToChromArrays2.MIN_SCAFFOLD+
				" merge="+FastaToChromArrays2.MERGE_SCAFFOLDS+" noReferenceTree="+noReferenceTree+" pass="+pass);
		return pass;
	}

	/** Four literal bases, exact retained bounds and literal in-memory metadata form the success oracle. */
	private static boolean tinyReference(ArrayList<ChromosomeArray> result){
		if(result==null || result.size()!=1){return false;}
		final ChromosomeArray ca=result.get(0);
		return ca.minIndex==0 && ca.maxIndex==3 && Arrays.equals(ca.array, MOTIF) &&
				FastaToChromArrays2.INFO_LIST!=null && FastaToChromArrays2.INFO_LIST.contains("1\t1\t1\t4\t4\t0\t0\t0") &&
				FastaToChromArrays2.SCAF_LIST!=null && FastaToChromArrays2.SCAF_LIST.contains("1\t1\t0\t4\tone");
	}

	/** Forces a metadata close failure while the independent scaffold writer succeeds.
	 * Setup failures occur outside the converter catch and cannot count as rejection. */
	private static void metadata(Path work) throws Exception{
		final byte[] raw=payload(), compressed=gzip(raw);
		verifyGzip(compressed, raw, false);
		final Path input=work.resolve("reference.fa.gz");
		Files.write(input, compressed);
		final Path genome=work.resolve("ref/genome/17"), full=Paths.get("/dev/full");
		require(Files.exists(full), "Metadata fault fixture requires /dev/full");
		Files.createDirectories(genome);
		final Path info=genome.resolve("info.txt"), summary=genome.resolve("summary.txt");
		Files.createSymbolicLink(info, full);
		require(!Files.exists(summary), "Fresh fixture must not contain a completion summary");
		final String[] convert={"in="+input, "build=17", "path="+work,
				"nodisk=f", "writechroms=f", "retain=t", "genscaffoldinfo=t",
				"startpad=0", "midpad=0", "stoppad=0", "minscaf=1", "mergescaffolds=t"};
		RuntimeException failure=null;
		try{FastaToChromArrays2.main2(convert);}
		catch(RuntimeException problem){failure=problem;}
		System.out.println("metadataFailure="+failure+" summaryExists="+Files.exists(summary));
		require(failure!=null, "Reference converter returned normally after metadata output failure");
		require(failure.getMessage()!=null && failure.getMessage().contains("reference metadata"),
				"Converter failed before the checked metadata-completion branch: "+failure);
		require(!Files.exists(summary), "Failed metadata output published a completion summary");
		require(Files.isSymbolicLink(info) && Files.readSymbolicLink(info).equals(full),
				"Converter replaced the output-failure target");
		final String scaffolds=ReadWrite.readString(genome.resolve("scaffolds.txt.gz").toString());
		require(scaffolds.contains("1\t1\t0\t100000\tone\n"),
				"Independent scaffold writer did not finish its literal record");
	}

	/** Fixed FASTA with short lines, one contig, no ambiguity or padding. */
	private static byte[] payload(){
		final ByteArrayOutputStream text=new ByteArrayOutputStream();
		for(byte b : ">one\n".getBytes(StandardCharsets.US_ASCII)){text.write(b);}
		for(int i=0; i<BASES; i++){
			text.write(MOTIF[i%MOTIF.length]);
			if(i%80==79){text.write('\n');}
		}
		return text.toByteArray();
	}

	private static byte[] gzip(byte[] raw) throws Exception{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(GZIPOutputStream out=new GZIPOutputStream(bytes)){out.write(raw);}
		return bytes.toByteArray();
	}

	/** Independently establishes correct payload bytes or a CRC-specific failure. */
	private static void verifyGzip(byte[] bytes, byte[] expected, boolean corrupt) throws Exception{
		final ByteArrayOutputStream decoded=new ByteArrayOutputStream();
		ZipException failure=null;
		try(GZIPInputStream in=new GZIPInputStream(new ByteArrayInputStream(bytes))){
			final byte[] buffer=new byte[8192];
			for(int n; (n=in.read(buffer))!=-1;){decoded.write(buffer, 0, n);}
		}catch(ZipException problem){failure=problem;}
		require(corrupt ? failure!=null : failure==null && Arrays.equals(decoded.toByteArray(), expected),
				"Independent gzip control failed");
	}

	private static void require(boolean ok, String message){
		if(!ok){throw new AssertionError(message);}
	}
	private static final int BASES=100000;
	private static final byte[] MOTIF={'A', 'C', 'G', 'T'};
}
