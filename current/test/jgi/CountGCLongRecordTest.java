package test.jgi;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;

import jgi.CountGC;
import parse.Parse;

/** Baseline/fixed streaming CountGC per-record bucket overflow controls. */
public final class CountGCLongRecordTest{

	public static void main(final String[] args) throws Exception{
		for(final String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase();
			final String value=split.length>1 ? split[1] : null;
			if("phase".equals(key)){
				phase=value;
			}else if("outdir".equals(key)){
				outdir=value;
			}else{
				throw new AssertionError("Unknown argument: "+arg);
			}
		}
		check("baseline".equals(phase) || "fixed".equals(phase),
				"phase must be baseline|fixed: "+phase);
		check(outdir!=null, "outdir is required");

		final Path scratch=Files.createDirectory(Paths.get(outdir).toAbsolutePath());
		final Case[] cases={
				new Case("fasta_short", true, SHORT_BASES, 21),
				new Case("fastq_short", false, SHORT_BASES, 33),
				new Case("fasta_long", true, LONG_BASES, 2147483660L),
				new Case("fastq_long", false, LONG_BASES, 4294967312L)
		};
		writeConfig(scratch, cases);
		writeExpected(scratch, cases);
		final StringBuilder observations=new StringBuilder("label\tstdout\tcounts\tbytes\tbases\tqualities\n");
		for(final Case testCase : cases){
			checkCase(scratch, observations, testCase);
		}
		Files.write(scratch.resolve("observed.tsv"), observations.toString().getBytes(StandardCharsets.UTF_8));
		System.out.println("PASS CountGCLongRecordTest phase="+phase+" checks="+checks);
	}

	private static void checkCase(final Path scratch, final StringBuilder observations,
			final Case testCase) throws Exception{
		final GeneratedLongRecordInputStream stream=new GeneratedLongRecordInputStream(
				testCase.label, testCase.fasta, testCase.bases);
		final Captured actual=capture(testCase, stream);
		final String expectedStdout=expectedStdout(testCase);
		final long[] expectedCounts=expectedCounts(testCase);
		final long expectedQualities=testCase.fasta ? 0 : testCase.bases;

		Files.write(scratch.resolve(testCase.label+".expected.stdout"),
				expectedStdout.getBytes(StandardCharsets.UTF_8));
		Files.write(scratch.resolve(testCase.label+".observed.stdout"),
				actual.stdout.getBytes(StandardCharsets.UTF_8));
		observations.append(testCase.label).append('\t')
				.append(escape(actual.stdout)).append('\t')
				.append(Arrays.toString(actual.counts)).append('\t')
				.append(stream.bytesEmitted()).append('\t')
				.append(stream.basesEmitted()).append('\t')
				.append(stream.qualitiesEmitted()).append('\n');
		Files.write(scratch.resolve("observed.tsv"),
				observations.toString().getBytes(StandardCharsets.UTF_8));

		check(expectedStdout.equals(actual.stdout), testCase.label+" stdout differed: "+actual.stdout);
		check(Arrays.equals(expectedCounts, actual.counts),
				testCase.label+" counts were "+Arrays.toString(actual.counts));
		check(stream.bytesEmitted()==testCase.expectedBytes,
				testCase.label+" bytes "+stream.bytesEmitted()+" != "+testCase.expectedBytes);
		check(stream.basesEmitted()==testCase.bases,
				testCase.label+" bases "+stream.basesEmitted()+" != "+testCase.bases);
		check(stream.qualitiesEmitted()==expectedQualities,
				testCase.label+" qualities "+stream.qualitiesEmitted()+" != "+expectedQualities);
	}

	private static Captured capture(final Case testCase,
			final GeneratedLongRecordInputStream stream) throws Exception{
		final int oldFormat=CountGC.FORMAT;
		final boolean oldSummary=CountGC.SUMMARY_ONLY;
		final PrintStream oldOut=System.out;
		final ByteArrayOutputStream capturedBytes=new ByteArrayOutputStream();
		final PrintStream captured=new PrintStream(capturedBytes, true,
				StandardCharsets.UTF_8.name());
		try{
			CountGC.FORMAT=1;
			CountGC.SUMMARY_ONLY=false;
			System.setOut(captured);
			final long[] counts=testCase.fasta ? CountGC.countFasta(stream, null) :
				CountGC.countFastq(stream, null);
			captured.flush();
			return new Captured(capturedBytes.toString(StandardCharsets.UTF_8.name()), counts);
		}finally{
			System.setOut(oldOut);
			CountGC.FORMAT=oldFormat;
			CountGC.SUMMARY_ONLY=oldSummary;
			captured.close();
		}
	}

	private static void writeConfig(final Path scratch, final Case[] cases) throws IOException{
		final StringBuilder sb=new StringBuilder("label\tformat\tbases\tqualities\texpectedBytes\n");
		for(final Case testCase : cases){
			sb.append(testCase.label).append('\t')
					.append(testCase.fasta ? "FASTA" : "FASTQ").append('\t')
					.append(testCase.bases).append('\t')
					.append(testCase.fasta ? 0 : testCase.bases).append('\t')
					.append(testCase.expectedBytes).append('\n');
		}
		Files.write(scratch.resolve("config.tsv"), sb.toString().getBytes(StandardCharsets.UTF_8));
	}

	private static void writeExpected(final Path scratch, final Case[] cases) throws IOException{
		final StringBuilder sb=new StringBuilder("label\tstdout\tcounts\tbytes\tbases\tqualities\n");
		for(final Case testCase : cases){
			sb.append(testCase.label).append('\t')
					.append(escape(expectedStdout(testCase))).append('\t')
					.append(Arrays.toString(expectedCounts(testCase))).append('\t')
					.append(testCase.expectedBytes).append('\t')
					.append(testCase.bases).append('\t')
					.append(testCase.fasta ? 0 : testCase.bases).append('\n');
		}
		Files.write(scratch.resolve("expected.tsv"), sb.toString().getBytes(StandardCharsets.UTF_8));
	}

	private static String expectedStdout(final Case testCase){
		final boolean wrap="baseline".equals(phase) && testCase.bases>Integer.MAX_VALUE;
		final long length=wrap ? Integer.MIN_VALUE : testCase.bases;
		final String a=wrap ? "-2147483648.00000" : "1.00000";
		final String row=testCase.label+'\t'+length+'\t'+a
				+"\t0.00000\t0.00000\t0.00000\t0.00000\t0.00000\n";
		return testCase.fasta ? row+'\n' : row;
	}

	private static long[] expectedCounts(final Case testCase){
		final long a="baseline".equals(phase) && testCase.bases>Integer.MAX_VALUE ?
				Integer.MIN_VALUE : testCase.bases;
		return new long[] {a, 0, 0, 0, 0, 0};
	}

	private static String escape(final String s){
		return s.replace("\\", "\\\\").replace("\t", "\\t").replace("\r", "\\r")
				.replace("\n", "\\n");
	}

	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
		checks++;
	}

	private static final class GeneratedLongRecordInputStream extends InputStream{
		GeneratedLongRecordInputStream(final String label, final boolean fasta_, final long bases_){
			if(bases_<0){throw new IllegalArgumentException("bases must be nonnegative");}
			basesRemaining=bases_;
			qualitiesRemaining=fasta_ ? 0 : bases_;
			prefix=((fasta_ ? ">" : "@")+label+'\n').getBytes(StandardCharsets.US_ASCII);
			separator=fasta_ ? EMPTY : "\n+\n".getBytes(StandardCharsets.US_ASCII);
			suffix=fasta_ ? EMPTY : "\n".getBytes(StandardCharsets.US_ASCII);
		}

		@Override
		public int read() throws IOException{
			final int len=read(single, 0, 1);
			return len<1 ? -1 : single[0]&0xFF;
		}

		@Override
		public int read(final byte[] b, final int off, final int len){
			if(b==null){throw new NullPointerException();}
			if(off<0 || len<0 || len>b.length-off){throw new IndexOutOfBoundsException();}
			if(len==0){return 0;}
			int written=0;
			while(written<len && stage<5){
				if(stage==0){written+=copyArray(prefix, b, off+written, len-written);}
				else if(stage==1){written+=fillRun((byte)'A', b, off+written, len-written, true);}
				else if(stage==2){written+=copyArray(separator, b, off+written, len-written);}
				else if(stage==3){written+=fillRun((byte)'I', b, off+written, len-written, false);}
				else if(stage==4){written+=copyArray(suffix, b, off+written, len-written);}
			}
			return written>0 ? written : -1;
		}

		private int copyArray(final byte[] src, final byte[] dest, final int off, final int len){
			if(arrayPos>=src.length){arrayPos=0; stage++; return 0;}
			final int toCopy=Math.min(len, src.length-arrayPos);
			System.arraycopy(src, arrayPos, dest, off, toCopy);
			arrayPos+=toCopy;
			bytes+=toCopy;
			return toCopy;
		}

		private int fillRun(final byte value, final byte[] dest, final int off,
				final int len, final boolean bases){
			long remaining=bases ? basesRemaining : qualitiesRemaining;
			if(remaining<1){stage++; return 0;}
			final int toCopy=(int)Math.min(remaining, len);
			Arrays.fill(dest, off, off+toCopy, value);
			if(bases){basesRemaining-=toCopy; basesEmitted+=toCopy;}
			else{qualitiesRemaining-=toCopy; qualitiesEmitted+=toCopy;}
			bytes+=toCopy;
			return toCopy;
		}

		long bytesEmitted(){return bytes;}
		long basesEmitted(){return basesEmitted;}
		long qualitiesEmitted(){return qualitiesEmitted;}

		private final byte[] prefix;
		private final byte[] separator;
		private final byte[] suffix;
		private long basesRemaining;
		private long qualitiesRemaining;
		private long bytes=0;
		private long basesEmitted=0;
		private long qualitiesEmitted=0;
		private int stage=0;
		private int arrayPos=0;
		private final byte[] single=new byte[1];
		private static final byte[] EMPTY=new byte[0];
	}

	private static final class Case{
		Case(final String label_, final boolean fasta_, final long bases_,
				final long expectedBytes_){
			label=label_;
			fasta=fasta_;
			bases=bases_;
			expectedBytes=expectedBytes_;
		}

		final String label;
		final boolean fasta;
		final long bases;
		final long expectedBytes;
	}

	private static final class Captured{
		Captured(final String stdout_, final long[] counts_){
			stdout=stdout_;
			counts=counts_;
		}

		final String stdout;
		final long[] counts;
	}

	private static String phase=null;
	private static String outdir=null;
	private static int checks=0;
	private static final long SHORT_BASES=8;
	private static final long LONG_BASES=1L+Integer.MAX_VALUE;
}
