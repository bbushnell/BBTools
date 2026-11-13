package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Locale;

import barcode.NovaDemux;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Finite NovaDemux regressions for mapin legacy counts, explicit FASTQ metadata,
 * and subset-restricted routing with optional legacy statistics.
 * Baseline mode requires the current exception paths for BAR008/BAR009; fixed mode
 * requires those exact inputs to demultiplex like the adjacent counted and
 * extensionless controls.
 * @author Brian Bushnell, Jean
 */
public final class NovaDemuxMapinMetadataTest{

	/** Requires mode=baseline|fixed and a supported case plus a new dir=path. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, NovaDemuxMapinMetadataTest.class, false);
		try{
			String mode=null, testCase=null, directory=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("case")){testCase=value;}
				else if(key.equals("dir")){directory=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			if(directory==null){throw new IllegalArgumentException("Set dir to a nonexistent output directory.");}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("NovaDemux mapin/metadata regression requires -ea.");}

			final Path dir=Files.createDirectory(Paths.get(directory));
			final boolean baseline="baseline".equals(mode);
			if("mapin".equals(testCase)){
				checkMapin(dir, baseline);
			}else if("counted".equals(testCase)){
				checkCountedLegacy(dir);
			}else if("metadata".equals(testCase)){
				checkMetadataFastq(dir, baseline);
			}else if("extensionless".equals(testCase)){
				checkExtensionlessMetadata(dir);
			}else if("subset-legacy-control".equals(testCase)){
				checkSubsetRouting(dir, baseline, false, false);
			}else if("subset-legacy".equals(testCase)){
				checkSubsetRouting(dir, baseline, true, false);
			}else if("subset-multiwriter-control".equals(testCase)){
				checkSubsetRouting(dir, baseline, false, true);
			}else if("subset-multiwriter".equals(testCase)){
				checkSubsetRouting(dir, baseline, true, true);
			}else{
				throw new IllegalArgumentException("Set case=mapin, case=counted, case=metadata, case=extensionless, "+
						"case=subset-legacy-control, case=subset-legacy, case=subset-multiwriter-control, "+
						"or case=subset-multiwriter.");
			}
			System.out.println("NOVA_MAPIN_METADATA mode="+mode+" case="+testCase+" failed=0");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; one fresh JVM owns each filesystem fixture. */
	private NovaDemuxMapinMetadataTest(){}

	/** Reproduces mapin+legacypath's null legacy count collection, or the proposed counted repair. */
	private static void checkMapin(final Path dir, final boolean baseline) throws IOException{
		writeTwoBarcodeInput(dir.resolve("reads"));
		writeMap(dir);
		Files.createDirectory(dir.resolve("legacy"));
		boolean npe=false;
		try{NovaDemux.main(mapinCommand(dir));}
		catch(final NullPointerException e){
			if(!baseline || !stackContains(e, "barcode.LegacyFileWriter", "writeTopUnknownBarcodes")){
				throw e;
			}
			npe=true;
		}
		check(npe==baseline, "mapin legacy NullPointerException: expected="+baseline+", observed="+npe);
		if(!baseline){
			waitForLegacyWriters(dir.resolve("legacy"));
			checkLegacyOutputs(dir);
		}
		System.out.println("mapin_legacy_npe="+npe);
	}

	/** Confirms ordinary input counting still supplies legacy reports on the same reads. */
	private static void checkCountedLegacy(final Path dir) throws IOException{
		writeTwoBarcodeInput(dir.resolve("reads"));
		Files.createDirectory(dir.resolve("legacy"));
		NovaDemux.main(countedCommand(dir));
		waitForLegacyWriters(dir.resolve("legacy"));
		checkLegacyOutputs(dir);
	}

	/** Reproduces suppressed .fq prefix probing after explicit interleaved=f, or the proposed metadata repair. */
	private static void checkMetadataFastq(final Path dir, final boolean baseline) throws IOException{
		final Path input=dir.resolve("reads.fq");
		final Path mapOut=dir.resolve("mapout.tsv");
		writeOneBarcodeInput(input);
		boolean assertion=false;
		try{NovaDemux.main(metadataCommand(input, mapOut));}
		catch(final AssertionError e){
			if(!baseline || !stackContains(e, "barcode.PCRMatrix", "<init>")){
				throw e;
			}
			assertion=true;
		}
		check(assertion==baseline, "explicit .fq metadata AssertionError: expected="+baseline+
				", observed="+assertion);
		if(baseline){check(!Files.exists(mapOut), "Metadata baseline must fail before writing mapout");}
		else{checkAssignmentMap(mapOut);}
		System.out.println("fastq_metadata_assertion="+assertion);
	}

	/** Confirms an extensionless FASTQ is still probed under the same explicit interleaving flag. */
	private static void checkExtensionlessMetadata(final Path dir) throws IOException{
		final Path input=dir.resolve("reads");
		final Path mapOut=dir.resolve("mapout.tsv");
		writeOneBarcodeInput(input);
		NovaDemux.main(metadataCommand(input, mapOut));
		checkAssignmentMap(mapOut);
	}

	/** Excluded assigned reads must not retain legacy-statistics destinations into matched output. */
	private static void checkSubsetRouting(final Path dir, final boolean baseline, final boolean legacy,
			final boolean multiWriter) throws IOException{
		writeTwoBarcodeInput(dir.resolve("reads.fq"));
		writeSubsetMap(dir);
		if(legacy){Files.createDirectory(dir.resolve("legacy"));}
		NovaDemux.main(subsetCommand(dir, legacy, multiWriter));
		if(legacy){waitForLegacyWriters(dir.resolve("legacy"));}
		checkText("Included subset output",
				"@INST:1:FC:1:1101:100:100 1:N:0:ACGT\n"+
				"ACGT\n+\nIIII\n", read(dir.resolve("out_ACGT.fq")));
		final Path excluded=dir.resolve("out_TGCA.fq");
		if(baseline && legacy){
			checkText("Leaked excluded subset output",
					"@INST:1:FC:1:1101:100:101 1:N:0:TGCA\n"+
					"TGCA\n+\nIIII\n", read(excluded));
		}else{
			check(!Files.exists(excluded),
					"Excluded assigned barcode must not retain a destination into matched output");
		}
		assertNoPayload(dir.resolve("unknown.fq"));
	}

	/** One assigned read and one unknown barcode for legacy map reports. */
	private static void writeTwoBarcodeInput(final Path path) throws IOException{
		Files.write(path, ("@INST:1:FC:1:1101:100:100 1:N:0:ACGT\n"+
				"ACGT\n+\nIIII\n"+
				"@INST:1:FC:1:1101:100:101 1:N:0:TGCA\n"+
				"TGCA\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** A single assigned read is sufficient for metadata-driven matrix construction. */
	private static void writeOneBarcodeInput(final Path path) throws IOException{
		Files.write(path, ("@INST:1:FC:1:1101:100:100 1:N:0:ACGT\n"+
				"ACGT\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** One explicit assignment map: the observed ACGT barcode maps to the expected ACGT sample. */
	private static void writeMap(final Path dir) throws IOException{
		Files.write(dir.resolve("map.tsv"), "ACGT\tACGT\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** Assignment map with one included and one explicitly excluded assigned destination. */
	private static void writeSubsetMap(final Path dir) throws IOException{
		Files.write(dir.resolve("map.tsv"), "ACGT\tACGT\nTGCA\tTGCA\n".getBytes(StandardCharsets.US_ASCII),
				StandardOpenOption.CREATE_NEW);
	}

	/** mapin takes the early assignment-map branch that used to skip legacy count setup. */
	private static String[] mapinCommand(final Path dir){
		return new String[] {
				"in="+dir.resolve("reads"),
				"expected=ACGT",
				"mapin="+dir.resolve("map.tsv"),
				"matrixtype=hdist",
				"legacypath="+dir.resolve("legacy"),
				"server=f",
				"interleaved=f",
				"qin=33",
				"ow=t",
				"t=1",
				"pigz=f",
				"unpigz=f",
				"loglog=f",
				"trimreadcomments=f",
		};
	}

	/** Ordinary counting is the nearby legacy-report control for the same reads and reports. */
	private static String[] countedCommand(final Path dir){
		return new String[] {
				"in="+dir.resolve("reads"),
				"expected=ACGT",
				"matrixtype=hdist",
				"legacypath="+dir.resolve("legacy"),
				"server=f",
				"interleaved=f",
				"qin=33",
				"ow=t",
				"t=1",
				"pigz=f",
				"unpigz=f",
				"loglog=f",
				"trimreadcomments=f",
		};
	}

	/** Explicit interleaved=f must not suppress barcode metadata needed to build the PCR matrix. */
	private static String[] metadataCommand(final Path input, final Path mapOut){
		return new String[] {
				"in="+input,
				"expected=ACGT",
				"mapout="+mapOut,
				"matrixtype=hdist",
				"server=f",
				"interleaved=f",
				"qin=33",
				"ow=t",
				"t=1",
				"pigz=f",
				"unpigz=f",
				"loglog=f",
				"trimreadcomments=f",
		};
	}

	/** Demultiplexes two assigned samples while outsubset permits only one of them to be written. */
	private static String[] subsetCommand(final Path dir, final boolean legacy, final boolean multiWriter){
		final ArrayList<String> args=new ArrayList<String>();
		args.add("in="+dir.resolve("reads.fq"));
		args.add("expected=ACGT,TGCA");
		args.add("outsubset=ACGT");
		args.add("mapin="+dir.resolve("map.tsv"));
		args.add("pattern="+dir.resolve("out_%.fq"));
		args.add("outu="+dir.resolve("unknown.fq"));
		args.add("matrixtype=hdist");
		if(legacy){args.add("legacypath="+dir.resolve("legacy"));}
		args.add("multiwriter="+multiWriter);
		args.add("writeempty=f");
		args.add("server=f");
		args.add("interleaved=f");
		args.add("qin=33");
		args.add("ow=t");
		args.add("t=1");
		args.add("pigz=f");
		args.add("unpigz=f");
		args.add("loglog=f");
		args.add("trimreadcomments=f");
		return args.toArray(new String[args.size()]);
	}

	/** Requires one exact headerless assignment row from PCRMatrix.printAssignmentMap. */
	private static void checkAssignmentMap(final Path mapOut) throws IOException{
		checkText("Assignment map", "ACGT\tACGT\t1\t0\n", read(mapOut));
	}

	/** Allows an absent or empty unknown file, but no rejected reads. */
	private static void assertNoPayload(final Path path) throws IOException{
		if(!Files.exists(path)){return;}
		checkText("Unknown subset output", "", read(path));
	}

	/** Verifies the single-index legacy reports reached after map construction. */
	private static void checkLegacyOutputs(final Path dir) throws IOException{
		final Path legacy=dir.resolve("legacy");
		checkText("Quality_Metrics",
				"Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30\n"+
				"0,Undetermined,ACGT,,1,4,4,160,40,1\n"+
				"0,Undetermined,,,1,4,4,160,40,1\n",
				read(legacy.resolve("Quality_Metrics.csv")));
		checkText("Demultiplex_Stats",
				"Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads\n"+
				"0,,ACGT,1,1,0,0,0.5000,1,0,0\n"+
				"0,Undetermined,,1,1,0,0,0.5000,1,0,0\n",
				read(legacy.resolve("Demultiplex_Stats.csv")));
		checkText("Top_Unknown_Barcodes",
				"Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads\n"+
				"0,TGCA,,1,1,0.500000\n",
				read(legacy.resolve("Top_Unknown_Barcodes.csv")));
		check(!Files.exists(legacy.resolve("Index_Hopping_Counts.csv")),
				"Single-index legacy output must not write Index_Hopping_Counts.csv");
	}

	/** Joins live legacy report writers left asynchronous by ByteStreamWriter.poison(). */
	private static void waitForLegacyWriters(final Path legacy){
		final String prefix="BSW for "+legacy.toString();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof ByteStreamWriter && thread.toString().startsWith(prefix)){
				final ByteStreamWriter bsw=(ByteStreamWriter)thread;
				bsw.waitForFinish();
				check(!bsw.errorState, "Legacy writer failed: "+bsw);
			}
		}
	}

	/** True iff the thrown call path contains the expected production frame. */
	private static boolean stackContains(final Throwable throwable, final String className, final String methodName){
		for(final StackTraceElement frame : throwable.getStackTrace()){
			if(frame.getClassName().equals(className) && frame.getMethodName().equals(methodName)){
				return true;
			}
		}
		return false;
	}

	/** Reads bounded fixture output as ASCII text. */
	private static String read(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	/** Requires byte-for-byte report text, including headers and newlines. */
	private static void checkText(final String label, final String expected, final String observed){
		if(!expected.equals(observed)){
			throw new AssertionError(label+" expected:\n"+expected+"observed:\n"+observed);
		}
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
}
