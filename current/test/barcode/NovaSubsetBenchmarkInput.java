package test.barcode;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Generates and independently verifies deterministic NovaDemux subset benchmark input and reports.
 * The input alternates the two Illumina headers from NovaDemuxMapinMetadataTest's subset fixture.
 * @author Jean
 */
public final class NovaSubsetBenchmarkInput{

	/** Requires mode=generate|verify|reports, path=FILE-or-DIR and an even reads count. */
	public static void main(final String[] args) throws IOException{
		final PreParser pp=new PreParser(args, NovaSubsetBenchmarkInput.class, false);
		try{
			String mode=null, path=null;
			long reads=-1;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(key.equals("path")){path=value;}
				else if(key.equals("reads")){reads=Long.parseLong(value);}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check(path!=null, "Set path to a FASTQ file for generate/verify or a legacy CSV directory for reports.");
			check(reads>0 && reads<=MAX_READS && (reads&1)==0,
					"Set an even reads count in [2,"+MAX_READS+"]; reads="+reads);
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Nova subset benchmark input requires -ea.");

			final long bytes;
			if("generate".equals(mode)){bytes=generate(Paths.get(path), reads);}
			else if("verify".equals(mode)){bytes=verify(Paths.get(path), reads);}
			else if("reports".equals(mode)){bytes=verifyReports(Paths.get(path), reads);}
			else{throw new IllegalArgumentException("Set mode=generate, mode=verify or mode=reports.");}
			System.out.println("NOVA_SUBSET_BENCHMARK_INPUT mode="+mode+" reads="+reads+
					" half="+(reads/2)+" tile=1101 bytes="+bytes);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; the launcher owns one deterministic file or legacy directory per invocation. */
	private NovaSubsetBenchmarkInput(){}

	/** Writes alternating ACGT/TGCA records with the exact subset-fixture headers and payloads. */
	private static long generate(final Path path, final long reads) throws IOException{
		check(!Files.exists(path), "Refusing to overwrite benchmark input: "+path);
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(), false, false, false);
		final ByteBuilder bb=new ByteBuilder(128);
		long bytes=0;
		writer.start();
		try{
			for(long i=0; i<reads; i++){
				final boolean acgt=((i&1)==0);
				bb.clear();
				appendRecord(bb, acgt);
				bytes+=bb.length();
				writer.print(bb);
			}
		}finally{
			if(writer.poisonAndWait()){throw new IOException("Benchmark input writer failed: "+path);}
		}
		check(Files.size(path)==bytes, "File size differs from the exact FASTQ byte count.");
		return bytes;
	}

	/** Verifies exact alternating bytes without line-ending normalization or FASTQ parsing lenience. */
	private static long verify(final Path path, final long reads) throws IOException{
		long bytes=0;
		try(final RawInput in=new RawInput(path)){
			for(long i=0; i<reads; i++){
				final boolean acgt=((i&1)==0);
				final byte[] expected=acgt ? ACGT_RECORD : TGCA_RECORD;
				for(final byte b : expected){in.expect(b, i);}
				bytes+=expected.length;
			}
			check(in.next()==-1, "Unexpected bytes after "+reads+" alternating Nova subset records.");
		}
		check(Files.size(path)==bytes, "File size differs from the exact FASTQ byte count.");
		return bytes;
	}

	/** Verifies four exact legacy CSV reports for an even half ACGT and half TGCA workload. */
	private static long verifyReports(final Path path, final long reads) throws IOException{
		final long half=reads/2, yield=4*half, quality=160*half;
		long bytes=0;
		bytes+=checkText(path.resolve("Demultiplex_Stats.csv"),
				"Lane,SampleID,Index,# Reads,# Perfect Index Reads,# One Mismatch Index Reads,# Two Mismatch Index Reads,% Reads,% Perfect Index Reads,% One Mismatch Index Reads,% Two Mismatch Index Reads\n"+
				"0,,ACGT,"+half+","+half+",0,0,0.5000,1,0,0\n"+
				"0,,TGCA,"+half+","+half+",0,0,0.5000,1,0,0\n"+
				"0,Undetermined,,0,0,0,0,0,0.0000,0.0000,0.0000\n");
		bytes+=checkText(path.resolve("Demultiplex_Tile_Stats.csv"),
				"Lane,SampleID,Index,Tile,# Reads,% Reads\n"+
				"0,null,ACGT,1101,"+half+",0.5000\n"+
				"0,null,TGCA,1101,"+half+",0.5000\n");
		bytes+=checkText(path.resolve("Quality_Metrics.csv"),
				"Lane,SampleID,index,index2,ReadNumber,Yield,YieldQ30,QualityScoreSum,Mean Quality Score (PF),% Q30\n"+
				"0,Undetermined,ACGT,,1,"+yield+","+yield+","+quality+",40,1\n"+
				"0,Undetermined,TGCA,,1,"+yield+","+yield+","+quality+",40,1\n"+
				"0,Undetermined,,,1,0,0,0,0.00,0.00\n");
		bytes+=checkText(path.resolve("Top_Unknown_Barcodes.csv"),
				"Lane,index,index2,# Reads,% of Unknown Barcodes,% of All Reads\n");
		return bytes;
	}

	/** Appends one of the two exact finite-fixture records. */
	private static void appendRecord(final ByteBuilder bb, final boolean acgt){
		bb.append(acgt ? ACGT_RECORD : TGCA_RECORD);
	}

	/** Reads a whole small legacy CSV and requires byte-for-byte expected text. */
	private static long checkText(final Path path, final String expected) throws IOException{
		final String observed=new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
		if(!expected.equals(observed)){
			throw new IOException("Unexpected legacy report "+path+"\nExpected:\n"+expected+"Observed:\n"+observed);
		}
		return observed.length();
	}

	/** Always enforces fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	/** Buffered raw-byte cursor: no CR stripping or implicit acceptance of a missing terminator. */
	private static final class RawInput implements AutoCloseable{
		RawInput(final Path path) throws IOException{input=Files.newInputStream(path, StandardOpenOption.READ);}
		int next() throws IOException{
			if(limit<0){return -1;}
			if(offset==limit){
				do{limit=input.read(buffer);}while(limit==0);
				offset=0;
				if(limit<0){return -1;}
			}
			return buffer[offset++]&255;
		}
		void expect(final byte expected, final long record) throws IOException{
			final int observed=next();
			if(observed!=(expected&255)){
				throw new IOException("Unexpected byte at record "+record+": expected="+(expected&255)+
						", observed="+observed);
			}
		}
		@Override
		public void close() throws IOException{input.close();}
		private final InputStream input;
		private final byte[] buffer=new byte[65536];
		private int offset=0, limit=0;
	}

	private static final long MAX_READS=32000000L;
	private static final byte[] ACGT_RECORD=("@INST:1:FC:1:1101:100:100 1:N:0:ACGT\n"+
			"ACGT\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII);
	private static final byte[] TGCA_RECORD=("@INST:1:FC:1:1101:100:101 1:N:0:TGCA\n"+
			"TGCA\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII);
}
