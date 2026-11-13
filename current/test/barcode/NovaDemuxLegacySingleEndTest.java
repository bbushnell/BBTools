package test.barcode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

import barcode.NovaDemux;
import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Single-end legacy-output regression for NovaDemux's absent-mate statistics path.
 * Baseline mode requires the current NullPointerException; fixed mode requires
 * one assigned read and one undetermined read to produce R1-only legacy reports
 * without empty R2 metric rows after asynchronous legacy writers finish.
 * The independently selected paired case preserves the existing R1/R2 quality rows.
 * @author Brian Bushnell, Jean
 */
public final class NovaDemuxLegacySingleEndTest{

	/** Requires mode=baseline|fixed, case=single|paired and a new dir=path; assertions must be enabled. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, NovaDemuxLegacySingleEndTest.class, false);
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
			if(!"single".equals(testCase) && !"paired".equals(testCase)){
				throw new IllegalArgumentException("Set case=single or case=paired.");
			}
			if(directory==null){throw new IllegalArgumentException("Set dir to a nonexistent output directory.");}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("NovaDemux legacy regression requires -ea.");}

			final Path dir=Files.createDirectory(Paths.get(directory));
			Files.createDirectory(dir.resolve("legacy"));
			if("single".equals(testCase)){
				writeInput(dir, "input");
				boolean npe=false;
				try{NovaDemux.main(singleEndCommand(dir));}
				catch(final NullPointerException e){
					if(!"baseline".equals(mode)){throw e;}
					boolean absentMate=false;
					for(final StackTraceElement frame : e.getStackTrace()){
						if(frame.getClassName().equals("barcode.LegacyFileWriter") && frame.getMethodName().equals("add")){
							absentMate=true;
							break;
						}
					}
					if(!absentMate){throw e;}
					npe=true;
				}
				if("baseline".equals(mode)){
					check(npe, "Baseline must fail on the absent single-end mate.");
				}else{
					check(!npe, "Fixed NovaDemux must skip the absent mate.");
					waitForLegacyWriters(dir.resolve("legacy"));
					checkSingleEndOutputs(dir);
				}
				System.out.println("NOVA_LEGACY_SINGLE mode="+mode+" case="+testCase+" npe="+npe);
			}else{
				writeInput(dir, "input1");
				writePairedMateInput(dir);
				NovaDemux.main(pairedCommand(dir));
				waitForLegacyWriters(dir.resolve("legacy"));
				checkPairedOutputs(dir);
				System.out.println("NOVA_LEGACY_SINGLE mode="+mode+" case="+testCase+" npe=false");
			}
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; one fresh JVM owns each filesystem fixture. */
	private NovaDemuxLegacySingleEndTest(){}

	/** Creates one assigned read and one undetermined read with Illumina-style barcode fields. */
	private static void writeInput(final Path dir, final String fname) throws IOException{
		Files.write(dir.resolve(fname), ("@INST:1:FC:1:1101:100:100 1:N:0:ACGT\n"+
				"ACGT\n+\nIIII\n"+
				"@INST:1:FC:1:1101:100:101 1:N:0:TGCA\n"+
				"TGCA\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Creates mates with distinct qualities so R2 legacy rows are verified independently. */
	private static void writePairedMateInput(final Path dir) throws IOException{
		Files.write(dir.resolve("input2"), ("@INST:1:FC:1:1101:100:100 2:N:0:ACGT\n"+
				"TGCA\n+\nHHHH\n"+
				"@INST:1:FC:1:1101:100:101 2:N:0:TGCA\n"+
				"ACGT\n+\nHHHH\n").getBytes(StandardCharsets.US_ASCII), StandardOpenOption.CREATE_NEW);
	}

	/** Routes one single-end expected barcode and one unexpected barcode. */
	private static String[] singleEndCommand(final Path dir){
		return new String[] {
				"in="+dir.resolve("input"),
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

	/** Routes one paired expected barcode and one unexpected barcode. */
	private static String[] pairedCommand(final Path dir){
		return new String[] {
				"in1="+dir.resolve("input1"),
				"in2="+dir.resolve("input2"),
				"expected=ACGT",
				"matrixtype=hdist",
				"legacypath="+dir.resolve("legacy"),
				"server=f",
				"qin=33",
				"ow=t",
				"t=1",
				"pigz=f",
				"unpigz=f",
				"loglog=f",
				"trimreadcomments=f",
		};
	}

	/** Verifies the R1-only legacy files that should exist for single-index input. */
	private static void checkSingleEndOutputs(final Path dir) throws IOException{
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

	/** Verifies paired legacy quality metrics retain their R1 and R2 rows. */
	private static void checkPairedOutputs(final Path dir) throws IOException{
		final Path legacy=dir.resolve("legacy");
		checkText("Paired Quality_Metrics",
				"Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30\n"+
				"0,Undetermined,ACGT,,1,4,4,160,40,1\n"+
				"0,Undetermined,ACGT,,2,4,4,156,39,1\n"+
				"0,Undetermined,,,1,4,4,160,40,1\n"+
				"0,Undetermined,,,2,4,4,156,39,1\n",
				read(legacy.resolve("Quality_Metrics.csv")));
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
