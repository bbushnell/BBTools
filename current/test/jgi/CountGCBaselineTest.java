package test.jgi;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import jgi.CountGC;
import parse.Parse;
import shared.Shared;

/** Finite CountGC baseline controls.
 * @author Jean, Sucrose
 */
public final class CountGCBaselineTest{

	public static void main(final String[] args) throws Exception{
		for(final String arg : args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String key=split[0].toLowerCase();
			final String value=split.length>1 ? split[1] : null;
			if("phase".equals(key)){
				phase=value;
			}else if("outdir".equals(key)){
				outdir=value;
			}else if("scope".equals(key)){
				scope=value;
			}else{
				throw new AssertionError("Unknown argument: "+arg);
			}
		}

		check(phase!=null, "phase is required");
		check(outdir!=null, "outdir is required");
		final boolean baseline="baseline".equals(phase);
		check(baseline || "fixed".equals(phase), "phase must be baseline or fixed: "+phase);
		if(!"default".equals(scope) && !"samepath".equals(scope)
				&& !"highbit".equals(scope) && !"writercleanup".equals(scope)
				&& !"fastqboundary".equals(scope) && !"outputerror".equals(scope)
				&& !"mainerrors".equals(scope) && !"accounting".equals(scope)
				&& !"inputownership".equals(scope) && !"readstringfqz".equals(scope)
				&& !"readobjectfqz".equals(scope)){
			throw new AssertionError("scope must be default, samepath, highbit, "
					+"writercleanup, fastqboundary, outputerror, mainerrors, "
					+"accounting, inputownership, readstringfqz, or readobjectfqz: "
					+scope);
		}

		final Path scratch=Files.createDirectory(Paths.get(outdir).toAbsolutePath());
		if("samepath".equals(scope)){
			checkSamePathInputOutput(scratch, baseline);
		}else if("highbit".equals(scope)){
			checkHighBitParsing(baseline);
		}else if("writercleanup".equals(scope)){
			checkWriterCleanup(scratch, baseline);
		}else if("fastqboundary".equals(scope)){
			checkSegmentedFastqCrLf(baseline);
		}else if("outputerror".equals(scope)){
			checkOutputErrorPropagation(scratch, baseline);
		}else if("mainerrors".equals(scope)){
			checkMainErrors(scratch, baseline);
		}else if("accounting".equals(scope)){
			checkAccounting(scratch, baseline);
		}else if("inputownership".equals(scope)){
			checkInputOwnership(scratch, baseline);
		}else if("readstringfqz".equals(scope)){
			checkPublicReadStringFqz(scratch);
		}else if("readobjectfqz".equals(scope)){
			checkPublicReadObjectFqz(scratch);
		}else{
			checkFastqInteriorAt(baseline);
			checkOutEqualsPath(scratch, baseline);
			checkIntGcNumeratorOverflow(baseline);
		}
		final String scopeString="default".equals(scope) ? "" : " scope="+scope;
		System.out.println("PASS CountGCBaselineTest phase="+phase+scopeString+" checks="+checks);
	}

	private static void checkFastqInteriorAt(final boolean baseline) throws Exception{
		final int oldFormat=CountGC.FORMAT;
		final boolean oldSummary=CountGC.SUMMARY_ONLY;
		final PrintStream oldOut=System.out;
		final ByteArrayOutputStream capturedBytes=new ByteArrayOutputStream();
		final PrintStream captured=new PrintStream(capturedBytes, true, "US-ASCII");
		try{
			CountGC.FORMAT=1;
			CountGC.SUMMARY_ONLY=false;
			final byte[] input="@read@A\nACGTN\n+\nIIIII\n".getBytes(StandardCharsets.US_ASCII);
			System.setOut(captured);
			final long[] counts=CountGC.countFastq(new ByteArrayInputStream(input), null);
			captured.flush();
			check(Arrays.equals(new long[]{1, 1, 1, 1, 1, 0}, counts),
					"FASTQ aggregate counts were "+Arrays.toString(counts));
		} finally{
			System.setOut(oldOut);
			CountGC.FORMAT=oldFormat;
			CountGC.SUMMARY_ONLY=oldSummary;
			captured.close();
		}

		final String expected=(baseline ? "readA" : "read@A")
				+"\t5\t0.25000\t0.25000\t0.25000\t0.25000\t0.20000\t0.50000\n";
		check(expected.equals(capturedBytes.toString("US-ASCII")),
				"FASTQ record name output differed: "+capturedBytes.toString("US-ASCII"));
	}

	private static void checkHighBitParsing(final boolean baseline) throws Exception{
		final String asciiRow="plain\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final long[] asciiFastaCounts={1, 1, 1, 1, 1, 1};
		final long[] asciiFastqCounts={1, 1, 1, 1, 1, 0};
		checkCapture("ASCII FASTA", asciiFastaCounts, asciiRow+"\n",
				new CountCall(){
					@Override
					public long[] count() throws Exception{
						final byte[] input=">plain\nACGTN\n".getBytes(StandardCharsets.US_ASCII);
						return CountGC.countFasta(new ByteArrayInputStream(input), null);
					}
				});
		checkCapture("ASCII FASTQ", asciiFastqCounts, asciiRow,
				new CountCall(){
					@Override
					public long[] count() throws Exception{
						final byte[] input="@plain\nACGTN\n+\nIIIII\n"
								.getBytes(StandardCharsets.US_ASCII);
						return CountGC.countFastq(new ByteArrayInputStream(input), null);
					}
				});

		final String utfFastaRow="caf\u00E9\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n\n";
		if(baseline){
			checkArrayIndexFailure("baseline UTF-8 FASTA header", -87,
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=">caf\u00E9\nACGTN\n"
									.getBytes(StandardCharsets.UTF_8);
							return CountGC.countFasta(new ByteArrayInputStream(input), null);
						}
					});
		}else{
			checkCapture("fixed UTF-8 FASTA header", asciiFastaCounts, utfFastaRow,
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=">caf\u00E9\nACGTN\n"
									.getBytes(StandardCharsets.UTF_8);
							return CountGC.countFasta(new ByteArrayInputStream(input), null);
						}
					});
		}

		final String utfFastqRow="read\u00E9\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		checkCapture("UTF-8 FASTQ header", asciiFastqCounts,
				baseline ? utfFastqRow.substring(0, 4)+utfFastqRow.substring(5)
						: utfFastqRow,
				new CountCall(){
					@Override
					public long[] count() throws Exception{
						final byte[] input="@read\u00E9\nACGTN\n+\nIIIII\n"
								.getBytes(StandardCharsets.UTF_8);
						return CountGC.countFastq(new ByteArrayInputStream(input), null);
					}
				});

		if(baseline){
			checkArrayIndexFailure("baseline high-bit FASTA sequence", -61,
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=new byte[]{'>', 'b', 'a', 'd', '\n',
									'A', 'C', (byte)0xC3, 'T', '\n'};
							return CountGC.countFasta(new ByteArrayInputStream(input), null);
						}
					});
		}else{
			checkClearInvalidSequence("fixed high-bit FASTA sequence",
					"Invalid FASTA sequence byte 0xC3",
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=new byte[]{'>', 'b', 'a', 'd', '\n',
									'A', 'C', (byte)0xC3, 'T', '\n'};
							return CountGC.countFasta(new ByteArrayInputStream(input), null);
						}
					});
		}

		if(baseline){
			checkCapture("baseline high-bit FASTQ sequence",
					new long[]{1, 1, 0, 0, 0, 0},
					"bad\t2\t0.50000\t0.50000\t0.00000\t0.00000\t0.00000\t0.50000\n"
					+"IIII\t0\t0.00000\t0.00000\t0.00000\t0.00000\t0.00000\t0.00000\n\n",
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=new byte[]{'@', 'b', 'a', 'd', '\n',
									'A', 'C', (byte)0xC3, 'T', '\n',
									'+', '\n', 'I', 'I', 'I', 'I', '\n'};
							return CountGC.countFastq(new ByteArrayInputStream(input), null);
						}
					});
		}else{
			checkClearInvalidSequence("fixed high-bit FASTQ sequence",
					"Invalid FASTQ sequence byte 0xC3",
					new CountCall(){
						@Override
						public long[] count() throws Exception{
							final byte[] input=new byte[]{'@', 'b', 'a', 'd', '\n',
									'A', 'C', (byte)0xC3, 'T', '\n',
									'+', '\n', 'I', 'I', 'I', 'I', '\n'};
							return CountGC.countFastq(new ByteArrayInputStream(input), null);
						}
					});
		}
	}


	private static void checkSegmentedFastqCrLf(final boolean baseline) throws Exception{
		final String normalRow="one\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final String normalTwoRows=normalRow
				+"two\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final String zeroOneRow="one\t0\t0.00000\t0.00000\t0.00000\t0.00000"
				+"\t0.00000\t0.00000\n";
		final String zeroTwoRow="two\t0\t0.00000\t0.00000\t0.00000\t0.00000"
				+"\t0.00000\t0.00000\n";
		final String zeroPlusRow="+\t0\t0.00000\t0.00000\t0.00000\t0.00000"
				+"\t0.00000\t0.00000\n";
		final long[] normalCounts={1, 1, 1, 1, 1, 0};
		final long[] normalTwoCounts={2, 2, 2, 2, 2, 0};
		final long[] zeroCounts={0, 0, 0, 0, 0, 0};
		final byte[] crlfInput="@one\r\nACGTN\r\n+\r\nIIIII\r\n"
				.getBytes(StandardCharsets.US_ASCII);
		final byte[] twoCrlfInput=("@one\r\nACGTN\r\n+\r\nIIIII\r\n"
				+"@two\r\nACGTN\r\n+\r\nIIIII\r\n")
				.getBytes(StandardCharsets.US_ASCII);
		final byte[] lfInput="@one\nACGTN\n+\nIIIII\n".getBytes(StandardCharsets.US_ASCII);

		checkCapture("segmented CRLF FASTQ", baseline ? zeroCounts : normalCounts,
				baseline ? zeroOneRow+zeroPlusRow : normalRow, new CountCall(){
			@Override
			public long[] count() throws Exception{
				return CountGC.countFastq(new ChunkedInputStream(crlfInput, 1), null);
			}
		});

		checkCapture("two-record segmented CRLF FASTQ", baseline ? zeroCounts : normalTwoCounts,
				baseline ? zeroOneRow+zeroPlusRow+zeroTwoRow+zeroPlusRow : normalTwoRows,
				new CountCall(){
			@Override
			public long[] count() throws Exception{
				return CountGC.countFastq(new ChunkedInputStream(twoCrlfInput, 1), null);
			}
		});

		checkCapture("segmented LF FASTQ", normalCounts, normalRow, new CountCall(){
			@Override
			public long[] count() throws Exception{
				return CountGC.countFastq(new ChunkedInputStream(lfInput, 1), null);
			}
		});

		checkCapture("unsegmented CRLF FASTQ", normalCounts, normalRow, new CountCall(){
			@Override
			public long[] count() throws Exception{
				return CountGC.countFastq(new ByteArrayInputStream(crlfInput), null);
			}
		});
	}



	private static void checkOutputErrorPropagation(final Path scratch,
			final boolean baseline) throws Exception{
		requireDevFull();
		final String oneRow="one\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";

		checkWriterNormal("normal FASTA output control",
				scratch.resolve("outputerror-normal.fa.tsv"), oneRow,
				new long[]{1, 1, 1, 1, 1, 1}, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=">one\nACGTN\n".getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFasta(new ByteArrayInputStream(input), out);
			}
		});

		checkOutputErrorDevice("normal FASTQ output error", baseline,
				new long[]{1, 1, 1, 1, 1, 0}, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input="@one\nACGTN\n+\nIIIII\n".getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFastq(new ByteArrayInputStream(input), out);
			}
		});

		checkOutputErrorPrimary("parser FASTA output error", baseline,
				"Invalid FASTA sequence byte 0xC3", null, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=new byte[]{'>', 'o', 'k', '\n', 'A', 'C', 'G', 'T', 'N', '\n',
						'>', 'b', 'a', 'd', '\n', 'A', 'C', (byte)0xC3, 'T', '\n'};
				return CountGC.countFasta(new ByteArrayInputStream(input), out);
			}
		});

		final IOException fastqFailure=new IOException("synthetic FASTQ read failure");
		checkOutputErrorPrimary("read FASTQ output error", baseline,
				null, fastqFailure, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=("@ok\nACGTN\n+\nIIIII\n"
						+"@bad\nAC").getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFastq(new ThrowingInputStream(input, fastqFailure), out);
			}
		});
	}

	/** Checks exact parser byte totals and same-invocation borrowed-stdin benchmark rates. */
	private static void checkAccounting(final Path scratch, final boolean baseline)
			throws Exception{
		final String row="acct\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		checkAccountingParser("empty FASTA", true, "", new long[6], "", baseline ? -1 : 0);
		checkAccountingParser("empty FASTQ", false, "", new long[6], "", baseline ? -1 : 0);
		checkAccountingParser("one-record FASTA", true, ">acct\nACGTN\n",
				new long[]{1, 1, 1, 1, 1, 1}, row+"\n", baseline ? 11 : 12);
		checkAccountingParser("one-record FASTQ", false, "@acct\nACGTN\n+\nIIIII\n",
				new long[]{1, 1, 1, 1, 1, 0}, row, baseline ? 19 : 20);
		checkAccountingBenchmark(scratch, baseline);
	}

	/** Keeps byte accounting separate from the existing FASTA/FASTQ count/output differences. */
	private static void checkAccountingParser(final String label, final boolean fasta,
			final String text, final long[] expectedCounts, final String expectedStdout,
			final long expectedBytes) throws Exception{
		final long oldLimsum=((Long)readStaticField(CountGC.class, "LIMSUM")).longValue();
		try{
			setStaticField(CountGC.class, "LIMSUM", Long.valueOf(0));
			final Captured actual=captureCount(new CountCall(){
				@Override
				public long[] count() throws Exception{
					final InputStream input=new ByteArrayInputStream(text.getBytes(StandardCharsets.US_ASCII));
					return fasta ? CountGC.countFasta(input, null) : CountGC.countFastq(input, null);
				}
			});
			final long bytes=((Long)readStaticField(CountGC.class, "LIMSUM")).longValue();
			recordObservation("accounting "+label, actual);
			System.out.println("OBS accounting "+label+"\tinputBytes="+text.length()+"\tLIMSUM="+bytes);
			check(actual.failure==null, label+" failed: "+exceptionSummary(actual.failure));
			check(Arrays.equals(expectedCounts, actual.counts),
					label+" counts differed: "+Arrays.toString(actual.counts));
			check(expectedStdout.equals(actual.stdout), label+" output differed: "+actual.stdout);
			check(bytes==expectedBytes, label+" byte total="+bytes+", expected="+expectedBytes);
		}finally{
			setStaticField(CountGC.class, "LIMSUM", Long.valueOf(oldLimsum));
		}
	}

	/** Compares formatted rates from one main call, without imposing a duration or cross-run rate. */
	private static void checkAccountingBenchmark(final Path scratch, final boolean baseline)
			throws Exception{
		check(!Files.exists(Paths.get("stdin.accounting.raw")),
				"stdin.accounting.raw must be absent so main borrows System.in");
		final char[] bases=new char[1<<20];
		Arrays.fill(bases, 'A');
		final TrackingInputStream input=new TrackingInputStream(new String(bases), null);
		final long oldLimsum=((Long)readStaticField(CountGC.class, "LIMSUM")).longValue();
		try{
			setStaticField(CountGC.class, "LIMSUM", Long.valueOf(0));
			final MainResult result=captureMain(input, "in=stdin.accounting.raw", "out=benchmark");
			saveMainObservation(scratch, "accounting-benchmark", result, input);
			recordMainObservation("accounting borrowed stdin benchmark", result, input);
			check(result.failure==null, "stdin benchmark failed: "+exceptionSummary(result.failure));
			check(result.stdout.isEmpty(), "stdin benchmark stdout differed: "+result.stdout);
			check(!input.closed, "stdin benchmark closed borrowed System.in");
			check(input.pos==(1<<20), "stdin benchmark consumed "+input.pos+" bytes instead of 1048576");
			final String raw=accountingSpeed(result.stderr, "Raw Speed:");
			final String uncompressed=accountingSpeed(result.stderr, "Uncompressed Speed:");
			final double rate=Double.parseDouble(uncompressed);
			check(rate>0 && !Double.isNaN(rate) && !Double.isInfinite(rate),
					"stdin benchmark lacks positive finite precision: "+uncompressed);
			check(baseline ? "0.00".equals(raw) : raw.equals(uncompressed),
					"stdin benchmark rates differed: raw="+raw+", uncompressed="+uncompressed);
			System.out.println("OBS accounting benchmark\tinputBytes=1048576\traw="+raw
					+"\tuncompressed="+uncompressed);
		}finally{
			setStaticField(CountGC.class, "LIMSUM", Long.valueOf(oldLimsum));
		}
	}

	/** Extracts exactly one printed rate with its units; malformed or missing observations fail loudly. */
	private static String accountingSpeed(final String stderr, final String prefix){
		String rate=null;
		for(final String line : stderr.split("\\r?\\n")){
			if(line.startsWith(prefix)){
				if(rate!=null){throw new AssertionError("Duplicate benchmark rate "+prefix+": "+stderr);}
				final String[] fields=line.substring(prefix.length()).trim().split("\\s+");
				if(fields.length!=2 || !"MBytes/s".equals(fields[1])){
					throw new AssertionError("Malformed benchmark rate "+prefix+": "+line);
				}
				rate=fields[0];
			}
		}
		if(rate==null){throw new AssertionError("Missing benchmark rate "+prefix+": "+stderr);}
		return rate;
	}

	private static void checkMainErrors(final Path scratch, final boolean baseline)
			throws Exception{
		check(!Files.exists(Paths.get("stdin.mainerrors.fa")),
				"stdin.mainerrors.fa fixture must resolve through raw System.in");
		checkMainCliOutputFailure(scratch, baseline);
		checkMainBorrowedStdinNormal(scratch);
		checkMainBorrowedStdinFailure(scratch, baseline);
		if(!baseline){checkCloseOwnedInputContracts();}
	}

	private static void checkMainCliOutputFailure(final Path scratch,
			final boolean baseline) throws Exception{
		requireDevFull();
		final Path input=scratch.resolve("mainout.fq");
		final Path stdout=scratch.resolve("mainerrors-output.stdout");
		final Path stderr=scratch.resolve("mainerrors-output.stderr");
		Files.write(input, "@one\nACGTN\n+\nIIIII\n".getBytes(StandardCharsets.US_ASCII));

		final Path countgc=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.exists(countgc), "countgc.sh not found at "+countgc);
		final int exit=run(scratch, "mainerrors-output", stdout, stderr, "bash",
				countgc.toString(), "-Xmx128m", "-Xms32m", "in="+input.toString(),
				"out="+DEV_FULL, "format=1", "ow=t", "simd=f");
		Files.write(scratch.resolve("mainerrors-output.exit"),
				Arrays.asList(Integer.toString(exit)), StandardCharsets.US_ASCII);

		final String stdoutText=readString(stdout);
		final String stderrText=readString(stderr);
		final String outputMessage="java.io.IOException: CountGC output error for "+DEV_FULL;
		final int outputIndex=stderrText.indexOf(outputMessage);
		final int npeIndex=stderrText.indexOf("java.lang.NullPointerException");
		check(exit!=0, "countgc.sh output-error exit was "+exit);
		check("".equals(stdoutText),
				"countgc.sh output-error wrote unexpected stdout: "+stdoutText);
		check(outputIndex>=0,
				"countgc.sh stderr lacked output IOException: "+stderrText);
		check(!stderrText.contains("Time:") && !stderrText.contains("Speed:"),
				"countgc.sh output-error reached timing output: "+stderrText);
		check(baseline ? npeIndex>outputIndex : npeIndex<0,
				"countgc.sh secondary NullPointerException mismatch: "+stderrText);
	}

	private static void checkMainBorrowedStdinNormal(final Path scratch)
			throws Exception{
		final Path output=scratch.resolve("mainerrors-stdin.tsv");
		final TrackingInputStream input=new TrackingInputStream(">stdin\nACGTN\n", null);
		final MainResult result=captureMain(input, "in=stdin.mainerrors.fa",
				"out="+output.toString(), "format=1", "ow=t");
		saveMainObservation(scratch, "mainerrors-stdin", result, input);
		recordMainObservation("main borrowed stdin normal", result, input);

		check(result.failure==null,
				"borrowed stdin main unexpectedly failed: "+exceptionSummary(result.failure));
		final String expectedRow="stdin\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final String actualOutput=readString(output);
		check(expectedRow.equals(actualOutput),
				"borrowed stdin output differed: "+actualOutput);
		check(!input.closed, "borrowed stdin was closed after normal main");
		check("".equals(result.stdout),
				"borrowed stdin main wrote unexpected stdout: "+result.stdout);
	}

	private static void checkMainBorrowedStdinFailure(final Path scratch,
			final boolean baseline)
			throws Exception{
		final IOException injected=new IOException("synthetic main stdin failure");
		final TrackingInputStream input=new TrackingInputStream(">ok\nACGTN\n>b", injected);
		final MainResult result=captureMain(input, "in=stdin.mainerrors.fa",
				"out=stdout", "format=1", "ow=t");
		saveMainObservation(scratch, "mainerrors-stdin-failure", result, input);
		recordMainObservation("main borrowed stdin failure", result, input);

		if(baseline){
			check(result.failure instanceof NullPointerException,
					"baseline escaping failure differed: "+exceptionSummary(result.failure));
		}else{
			check(result.failure==injected,
					"fixed escaping failure identity differed: "
					+exceptionSummary(result.failure));
		}
		check(baseline==result.stderr.contains("synthetic main stdin failure"),
				"stdin failure stderr mismatch: "+result.stderr);
		check(!input.closed, "borrowed stdin was closed after failing main");
		final String expectedRow="ok\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		check(expectedRow.equals(result.stdout),
				"failing borrowed stdin stdout differed: "+result.stdout);
	}

	private static void saveMainObservation(final Path scratch, final String stem,
			final MainResult result, final TrackingInputStream input) throws IOException{
		Files.write(scratch.resolve(stem+".stdout"),
				result.stdout.getBytes(StandardCharsets.UTF_8));
		Files.write(scratch.resolve(stem+".stderr"),
				result.stderr.getBytes(StandardCharsets.UTF_8));
		Files.write(scratch.resolve(stem+".failure"),
				Arrays.asList(exceptionSummary(result.failure)), StandardCharsets.US_ASCII);
		Files.write(scratch.resolve(stem+".inputClosed"),
				Arrays.asList(Boolean.toString(input.closed)), StandardCharsets.US_ASCII);
		Files.write(scratch.resolve(stem+".inputPos"),
				Arrays.asList(Integer.toString(input.pos)), StandardCharsets.US_ASCII);
		Files.write(scratch.resolve(stem+".snapshot"),
				Arrays.asList(String.valueOf(result.snapshot)), StandardCharsets.US_ASCII);
	}

	private static void recordMainObservation(final String label,
			final MainResult result, final TrackingInputStream input){
		System.out.println("OBS "+label
				+"\tstdout="+escape(result.stdout)
				+"\tstderr="+escape(result.stderr)
				+"\texception="+exceptionSummary(result.failure)
				+"\tinputClosed="+input.closed);
	}

	private static void checkInputOwnership(final Path scratch, final boolean baseline)
			throws Exception{
		checkJavaGzipBorrowedStdin(scratch, baseline);
		checkExternalPigzCompletion(scratch, baseline);
		if(!baseline){
			checkManagedInputLifecycleContracts(scratch);
			checkPublicFqzWholeFileRoutes(scratch);
		}
	}

	/** Proves compressed stdin aliases wrap borrowed System.in in the public Java gzip path. */
	private static void checkJavaGzipBorrowedStdin(final Path scratch,
			final boolean baseline) throws Exception{
		check(!Files.exists(Paths.get("stdin.java.fa.gz")),
				"stdin.java.fa.gz must be absent so main borrows compressed System.in");
		final Path output=scratch.resolve("inputownership-java-gzip.tsv");
		final byte[] compressed=gzipBytes(">java-gzip\nACGTN\n"
				.getBytes(StandardCharsets.US_ASCII));
		Files.write(scratch.resolve("inputownership-java-gzip.bytes.gz"), compressed);
		final TrackingInputStream input=new TrackingInputStream(compressed, null);
		final MainResult result=captureMain(input, "in=stdin.java.fa.gz",
				"out="+output.toString(), "format=1", "ow=t", "threads=1",
				"unpigz=f", "unbgzip=f", "gunzip=f", "nativebgzip=f");
		saveMainObservation(scratch, "inputownership-java-gzip", result, input);
		recordMainObservation("input ownership Java gzip stdin", result, input);

		check(result.failure==null,
				"Java gzip stdin main failed: "+exceptionSummary(result.failure));
		final String expectedRow="java-gzip\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final String actualOutput=readString(output);
		check(expectedRow.equals(actualOutput), "Java gzip output differed: "+actualOutput);
		check(input.pos==compressed.length, "Java gzip consumed "+input.pos
				+" bytes instead of "+compressed.length);
		check(input.closed==baseline, "Java gzip borrowed stdin closed="+input.closed
				+", baseline="+baseline);
		check("".equals(result.stdout),
				"Java gzip main wrote unexpected stdout: "+result.stdout);
		check(result.snapshot!=null && result.snapshot.isJavaGzipState(),
				"Java gzip backend flags differed: "+result.snapshot);
	}

	/** Proves the current public pigz route can ignore a decoder process failure. */
	private static void checkExternalPigzCompletion(final Path scratch,
			final boolean baseline) throws Exception{
		final Path launcher=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.isRegularFile(launcher), "countgc.sh not found at "+launcher);

		final Path evidence=Files.createDirectory(scratch.resolve("inputownership-pigz"));
		final Path evidenceBin=Files.createDirectory(evidence.resolve("bin"));
		final Path marker=evidence.resolve("pigz.marker");
		final Path fakePigz=evidenceBin.resolve("pigz");
		final Path input=evidence.resolve("input.fa.gz");
		final Path output=evidence.resolve("countgc.out");
		final Path stdout=evidence.resolve("countgc.stdout");
		final Path stderr=evidence.resolve("countgc.stderr");
		Files.write(input, gzipBytes(">ignored\nNNNN\n".getBytes(StandardCharsets.US_ASCII)));
		writeFakePigz(fakePigz, marker);

		final int exit=runWithPathFirst(scratch, "inputownership-pigz", stdout, stderr,
				evidenceBin, marker, "timeout", "--kill-after=5s", "30s",
				"bash", launcher.toString(), "-Xmx128m", "-Xms32m", "in="+input.toString(),
				"out="+output.toString(), "format=1", "ow=t", "unbgzip=f",
				"nativebgzip=f", "threads=2");
		Files.write(evidence.resolve("countgc.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);

		final String expectedMarker="pigz\t--version\n"
				+"pigz\t-c\t-d\t"+input.toString()+"\n";
		check(expectedMarker.equals(readString(marker)),
				"fake pigz argv marker differed: "+readString(marker));
		check("".equals(readString(stdout)),
				"pigz countgc.sh wrote unexpected stdout: "+readString(stdout));
		final String stderrText=readString(stderr);
		check(stderrText.contains("jgi.CountGC") && stderrText.contains("-Xmx128m")
				&& stderrText.contains("in="+input.toString()),
				"countgc.sh stderr lacked echoed Java command: "+stderrText);

		final String expectedRow="pigz\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		if(baseline){
			check(exit==0, "baseline pigz exit="+exit
					+"\nstdout:\n"+readString(stdout)+"\nstderr:\n"+stderrText);
			check(expectedRow.equals(readString(output)),
					"baseline pigz output differed: "+readString(output));
			check(stderrText.contains("Overall\t5\t0.25000\t0.25000\t0.25000\t0.25000"
					+"\t0.20000\t0.50000"),
					"baseline pigz stderr lacked Overall row: "+stderrText);
		}else{
			check(exit==1, "fixed pigz decoder failure exit was "+exit);
			check(expectedRow.equals(readString(output)),
					"fixed pigz output differed: "+readString(output));
			check(!stderrText.contains("Overall\t5\t0.25000\t0.25000\t0.25000\t0.25000"),
					"fixed pigz stderr reached Overall row: "+stderrText);
			check(!stderrText.contains("\nTime: \t"),
					"fixed pigz stderr reached Time report: "+stderrText);
			check(stderrText.contains("java.io.IOException: pigz decoder failed for "
					+input.toString()),
					"fixed pigz stderr lacked decoder failure: "+stderrText);
		}
	}

	/** Exercises private close-owned InputStream lifecycle without external tools. */
	private static void checkManagedInputLifecycleContracts(final Path scratch)
			throws Exception{
		final Class<?> helperClass=Class.forName("fileIO.ReadWrite$HelperInput");
		final Class<?> managedClass=Class.forName("fileIO.ReadWrite$ManagedInputStream");
		final Constructor<?> helperCtor=helperClass.getDeclaredConstructor();
		helperCtor.setAccessible(true);
		final Method closeNormally=helperClass.getDeclaredMethod("closeNormally");
		closeNormally.setAccessible(true);
		final Method closeAfterFailure=helperClass.getDeclaredMethod("closeAfterFailure",
				Throwable.class);
		closeAfterFailure.setAccessible(true);
		final Constructor<?> managedCtor=managedClass.getDeclaredConstructor(
				InputStream.class, helperClass);
		managedCtor.setAccessible(true);

		final byte[] payload="bulk".getBytes(StandardCharsets.US_ASCII);
		final BulkTrackingInputStream arrayInput=new BulkTrackingInputStream(payload);
		final InputStream bulkArray=managedInput(helperCtor, managedCtor, helperClass,
				arrayInput, null, null);
		final byte[] arrayBuffer=new byte[4];
		check(bulkArray.read(arrayBuffer)==4, "managed read(byte[]) length differed");
		check(Arrays.equals(payload, arrayBuffer),
				"managed read(byte[]) content differed: "+Arrays.toString(arrayBuffer));
		final BulkTrackingInputStream sliceInput=new BulkTrackingInputStream(payload);
		final InputStream bulkSlice=managedInput(helperCtor, managedCtor, helperClass,
				sliceInput, null, null);
		final byte[] sliceBuffer=new byte[6];
		check(bulkSlice.read(sliceBuffer, 2, 4)==4,
				"managed read(byte[], int, int) length differed");
		check(Arrays.equals(new byte[]{0, 0, 'b', 'u', 'l', 'k'}, sliceBuffer),
				"managed read(byte[], int, int) content differed: "
				+Arrays.toString(sliceBuffer));
		check(arrayInput.arrayReads==1 && arrayInput.sliceReads==0
				&& arrayInput.singleReads==0 && sliceInput.arrayReads==0
				&& sliceInput.sliceReads==1 && sliceInput.singleReads==0,
				"managed bulk reads were not directly delegated: array="
				+arrayInput.arrayReads+"/"+sliceInput.arrayReads+" slice="
				+arrayInput.sliceReads+"/"+sliceInput.sliceReads+" single="
				+arrayInput.singleReads+"/"+sliceInput.singleReads);

		final CloseTrackingInputStream normalClose=new CloseTrackingInputStream(null);
		final InputStream normal=managedInput(helperCtor, managedCtor, helperClass,
				normalClose, null, null);
		normal.close();
		normal.close();
		check(normalClose.closes==1, "managed close was not one-shot: "
				+normalClose.closes);

		final IOException closeFailure=new IOException("synthetic managed close failure");
		final CloseTrackingInputStream failingClose=
				new CloseTrackingInputStream(closeFailure);
		final String decoderName=scratch.resolve("managed-close.gz").toString();
		final InputStream failing=managedInput(helperCtor, managedCtor, helperClass,
				failingClose, decoderName, "pigz");
		addFinishedProcess(decoderName, 7);
		try{
			failing.close();
			throw new AssertionError("managed close failure was not thrown");
		}catch(IOException e){
			check(e==closeFailure,
					"managed close did not preserve stream failure identity: "
					+exceptionSummary(e));
			final Throwable[] suppressed=e.getSuppressed();
			check(suppressed.length==1
					&& suppressed[0] instanceof IOException
					&& suppressed[0].getMessage()!=null
					&& suppressed[0].getMessage().contains(
							"pigz decoder failed for "+decoderName),
					"managed decoder failure was not suppressed under close failure: "
					+Arrays.toString(suppressed));
		}
		failing.close();
		check(failingClose.closes==1, "managed failing close was not one-shot: "
				+failingClose.closes);

		final CloseTrackingInputStream readStringClose=
				new CloseTrackingInputStream(null);
		final String readStringName=scratch.resolve("readString-helper.fqz").toString();
		final Object readStringHelper=helperInput(helperCtor, helperClass,
				readStringClose, readStringClose, readStringName, "FQZ");
		addFinishedProcess(readStringName, 0);
		check(invokeHelperCloseNormally(closeNormally, readStringHelper)==null,
				"private readString-style helper close reported failure");
		check(readStringClose.closes==1,
				"private readString-style helper did not close exactly once: "
				+readStringClose.closes);

		final CloseTrackingInputStream readObjectClose=
				new CloseTrackingInputStream(null);
		final String readObjectName=scratch.resolve("readObject-helper.fqz").toString();
		final Object readObjectHelper=helperInput(helperCtor, helperClass,
				readObjectClose, readObjectClose, readObjectName, "FQZ");
		addFinishedProcess(readObjectName, 7);
		final Throwable readObjectFailure=
				invokeHelperCloseNormally(closeNormally, readObjectHelper);
		check(readObjectFailure instanceof IOException
				&& readObjectFailure.getMessage()!=null
				&& readObjectFailure.getMessage().contains(
						"FQZ decoder failed for "+readObjectName),
				"private readObject-style helper did not report FQZ completion: "
				+exceptionSummary(readObjectFailure));

		final IOException primary=new IOException("synthetic readString failure");
		final CloseTrackingInputStream readFailureClose=
				new CloseTrackingInputStream(null);
		final String readFailureName=
				scratch.resolve("readString-failure-helper.fqz").toString();
		final Object readFailureHelper=helperInput(helperCtor, helperClass,
				readFailureClose, readFailureClose, readFailureName, "FQZ");
		addFinishedProcess(readFailureName, 7);
		closeAfterFailure.invoke(readFailureHelper, primary);
		check(readFailureClose.closes==1
				&& primary.getSuppressed().length==1
				&& primary.getSuppressed()[0] instanceof IOException
				&& primary.getSuppressed()[0].getMessage()!=null
				&& primary.getSuppressed()[0].getMessage().contains(
						"FQZ decoder failed for "+readFailureName),
				"private read failure helper did not suppress FQZ completion: closes="
				+readFailureClose.closes+" suppressed="+exceptionSummaries(primary));
	}

	private static void checkPublicFqzWholeFileRoutes(final Path scratch)
			throws Exception{
		final Path evidence=Files.createDirectory(scratch.resolve("public-fqz"));
		checkPublicFqzReadString(evidence);
		checkPublicFqzReadObject(evidence);
	}

	private static void checkPublicFqzReadString(final Path evidence)
			throws Exception{
		final Path caseDir=Files.createDirectory(evidence.resolve("readString"));
		final Path bin=Files.createDirectory(caseDir.resolve("bin"));
		final Path marker=caseDir.resolve("fqz.marker");
		final Path fakeFqz=bin.resolve("fqz_comp");
		final Path child=caseDir.resolve("child");
		final Path stdout=caseDir.resolve("child.stdout");
		final Path stderr=caseDir.resolve("child.stderr");
		writeFakeWholeFileFqz(fakeFqz);

		final int exit=runWithPathFirst(evidence, "readString-fqz", stdout, stderr,
				bin, marker, "timeout", "--kill-after=5s", "30s",
				"bash", "dev/countgcbaselinetest.sh", "phase=fixed",
				"outdir="+child.toString(), "scope=readstringfqz");
		check(exit==0, "public readString FQZ child exit was "+exit
				+"\nstdout:\n"+readStringIfExists(stdout)
				+"\nstderr:\n"+readStringIfExists(stderr));
		final String expectedMarker="fqz\t-d\t\t"
				+child.resolve("normal-string.fqz").toString()+"\n"
				+"fqz\t-d\t\t"+child.resolve("exit7-string.fqz").toString()+"\n";
		check(expectedMarker.equals(readString(marker)),
				"public readString FQZ marker differed: "+readString(marker));
	}

	private static void checkPublicFqzReadObject(final Path evidence)
			throws Exception{
		final Path caseDir=Files.createDirectory(evidence.resolve("readObject"));
		final Path bin=Files.createDirectory(caseDir.resolve("bin"));
		final Path marker=caseDir.resolve("fqz.marker");
		final Path fakeFqz=bin.resolve("fqz_comp");
		final Path child=caseDir.resolve("child");
		final Path stdout=caseDir.resolve("child.stdout");
		final Path stderr=caseDir.resolve("child.stderr");
		writeFakeWholeFileFqz(fakeFqz);

		final int exit=runWithPathFirst(evidence, "readObject-fqz", stdout, stderr,
				bin, marker, "timeout", "--kill-after=5s", "30s",
				"bash", "dev/countgcbaselinetest.sh", "phase=fixed",
				"outdir="+child.toString(), "scope=readobjectfqz");
		check(exit==0, "public readObject FQZ child exit was "+exit
				+"\nstdout:\n"+readStringIfExists(stdout)
				+"\nstderr:\n"+readStringIfExists(stderr));
		final String expectedMarker="fqz\t-d\t\t"
				+child.resolve("normal-object.fqz").toString()+"\n"
				+"fqz\t-d\t\t"+child.resolve("exit7-object.fqz").toString()+"\n"
				+"fqz\t-d\t\t"+child.resolve("corrupt-object.fqz").toString()+"\n";
		check(expectedMarker.equals(readString(marker)),
				"public readObject FQZ marker differed: "+readString(marker));
	}

	private static void checkPublicReadStringFqz(final Path scratch){
		try{
			touch(scratch.resolve("normal-string.fqz"));
			touch(scratch.resolve("exit7-string.fqz"));
		}catch(IOException e){
			throw new AssertionError(e);
		}
		final String normal=scratch.resolve("normal-string.fqz").toString();
		final String exit7=scratch.resolve("exit7-string.fqz").toString();
		check("alpha\nbeta\n".equals(ReadWrite.readString(normal)),
				"public FQZ readString normal output differed");
		final RuntimeException failure=expectRuntime(new ThrowingCall(){
			@Override
			public void run(){
				ReadWrite.readString(exit7);
			}
		}, "public FQZ readString accepted exit7");
		check(causeMessageContains(failure, "FQZ decoder failed for "+exit7),
				"public FQZ readString did not report decoder failure: "
				+exceptionSummary(failure));
	}

	private static void checkPublicReadObjectFqz(final Path scratch){
		try{
			writeSerializedObject(scratch.resolve("object.payload"),
					"FQZ serialized payload");
			touch(scratch.resolve("normal-object.fqz"));
			touch(scratch.resolve("exit7-object.fqz"));
			touch(scratch.resolve("corrupt-object.fqz"));
		}catch(IOException e){
			throw new AssertionError(e);
		}
		final String normal=scratch.resolve("normal-object.fqz").toString();
		final String exit7=scratch.resolve("exit7-object.fqz").toString();
		final String corrupt=scratch.resolve("corrupt-object.fqz").toString();
		final Object object=ReadWrite.readObject(normal, true);
		check("FQZ serialized payload".equals(object),
				"public FQZ readObject normal object differed: "+object);
		final RuntimeException exit7Failure=expectRuntime(new ThrowingCall(){
			@Override
			public void run(){
				ReadWrite.readObject(exit7, true);
			}
		}, "public FQZ readObject accepted exit7");
		check(causeMessageContains(exit7Failure, "FQZ decoder failed for "+exit7),
				"public FQZ readObject did not report decoder failure: "
				+exceptionSummary(exit7Failure));

		final RuntimeException corruptFailure=expectRuntime(new ThrowingCall(){
			@Override
			public void run(){
				ReadWrite.readObject(corrupt, true);
			}
		}, "public FQZ readObject accepted a corrupt object stream");
		final Throwable cause=corruptFailure.getCause();
		check(cause instanceof IOException
				&& cause.getClass().getName().endsWith(".StreamCorruptedException"),
				"public FQZ readObject primary failure changed: "
				+exceptionSummary(corruptFailure));
		check(hasSuppressedMessage(cause, "FQZ decoder failed for "+corrupt),
				"public FQZ readObject corrupt stream did not suppress decoder "
				+"failure: "+exceptionSummaries(cause));
	}

	private static InputStream managedInput(final Constructor<?> helperCtor,
			final Constructor<?> managedCtor, final Class<?> helperClass,
			final InputStream inner, final String processName,
			final String processCodec) throws Exception{
		final Object helper=helperInput(helperCtor, helperClass, inner, inner,
				processName, processCodec);
		return (InputStream)managedCtor.newInstance(inner, helper);
	}

	private static Object helperInput(final Constructor<?> helperCtor,
			final Class<?> helperClass, final InputStream stream,
			final Closeable outer, final String processName,
			final String processCodec) throws Exception{
		final Object helper=helperCtor.newInstance();
		setField(helperClass, helper, "stream", stream);
		setField(helperClass, helper, "outer", outer);
		if(processName!=null){
			setField(helperClass, helper, "managed", Boolean.TRUE);
			setField(helperClass, helper, "processName", processName);
			setField(helperClass, helper, "processCodec", processCodec);
		}
		return helper;
	}

	private static Throwable invokeHelperCloseNormally(final Method closeNormally,
			final Object helper) throws Exception{
		try{
			closeNormally.invoke(helper);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	private static void addFinishedProcess(final String name, final int exit)
			throws Exception{
		final Method addProcess=ReadWrite.class.getDeclaredMethod("addProcess",
				String.class, Process.class);
		addProcess.setAccessible(true);
		addProcess.invoke(null, name, new FinishedProcess(exit));
	}

	private static void setField(final Class<?> owner, final Object instance,
			final String name, final Object value) throws Exception{
		final Field field=owner.getDeclaredField(name);
		field.setAccessible(true);
		field.set(instance, value);
	}

	private static void writeFakePigz(final Path fakePigz, final Path marker)
			throws IOException{
		final String script="#!/usr/bin/env bash\n"
				+"set -euo pipefail\n"
				+"marker=${FAKE_PIGZ_MARKER:?}\n"
				+"printf 'pigz' >> \"$marker\"\n"
				+"for arg in \"$@\"; do printf '\\t%s' \"$arg\" >> \"$marker\"; done\n"
				+"printf '\\n' >> \"$marker\"\n"
				+"if [[ $# -eq 1 && $1 = --version ]]; then\n"
				+"  printf 'pigz 2.3.1\\n'\n"
				+"  exit 0\n"
				+"fi\n"
				+"if [[ $# -eq 3 && $1 = -c && $2 = -d && -f $3 ]]; then\n"
				+"  printf '>pigz\\nACGTN\\n'\n"
				+"  exit 7\n"
				+"fi\n"
				+"printf 'unexpected pigz argv' >&2\n"
				+"exit 64\n";
		Files.write(fakePigz, script.getBytes(StandardCharsets.US_ASCII));
		check(fakePigz.toFile().setExecutable(true, false),
				"Could not make fake pigz executable at "+fakePigz);
	}

	private static void writeFakeWholeFileFqz(final Path fakeFqz) throws IOException{
		final String script="#!/usr/bin/env bash\n"
				+"set -euo pipefail\n"
				+"marker=${FAKE_PIGZ_MARKER:?}\n"
				+"printf 'fqz' >> \"$marker\"\n"
				+"for arg in \"$@\"; do printf '\\t%s' \"$arg\" >> \"$marker\"; done\n"
				+"printf '\\n' >> \"$marker\"\n"
				+"if [[ $# -eq 3 && $1 = -d && $2 = \"\" && -f $3 ]]; then\n"
				+"  base=${3##*/}\n"
				+"  case \"$base\" in\n"
				+"    normal-string.fqz)\n"
				+"      printf 'alpha\\nbeta\\n'\n"
				+"      exit 0\n"
				+"      ;;\n"
				+"    exit7-string.fqz)\n"
				+"      printf 'exit7\\n'\n"
				+"      exit 7\n"
				+"      ;;\n"
				+"    normal-object.fqz)\n"
				+"      cat \"$(dirname \"$3\")/object.payload\"\n"
				+"      exit 0\n"
				+"      ;;\n"
				+"    exit7-object.fqz)\n"
				+"      cat \"$(dirname \"$3\")/object.payload\"\n"
				+"      exit 7\n"
				+"      ;;\n"
				+"    corrupt-object.fqz)\n"
				+"      printf 'not a serialized object'\n"
				+"      exit 7\n"
				+"      ;;\n"
				+"  esac\n"
				+"fi\n"
				+"printf 'unexpected fqz_comp argv' >&2\n"
				+"exit 64\n";
		Files.write(fakeFqz, script.getBytes(StandardCharsets.US_ASCII));
		check(fakeFqz.toFile().setExecutable(true, false),
				"Could not make fake fqz_comp executable at "+fakeFqz);
	}

	private static void touch(final Path path) throws IOException{
		Files.write(path, new byte[0]);
	}

	private static void writeSerializedObject(final Path path, final Object object)
			throws IOException{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(ObjectOutputStream out=new ObjectOutputStream(bytes)){
			out.writeObject(object);
		}
		Files.write(path, bytes.toByteArray());
	}

	private static byte[] gzipBytes(final byte[] plain) throws IOException{
		final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
		try(final OutputStream gzip=new GZIPOutputStream(bytes)){
			gzip.write(plain);
		}
		return bytes.toByteArray();
	}

	private static void checkCloseOwnedInputContracts() throws Exception{
		final Method method=CountGC.class.getDeclaredMethod("closeOwnedInput",
				InputStream.class, Throwable.class);
		method.setAccessible(true);

		final CloseTrackingInputStream normal=new CloseTrackingInputStream(null);
		final Throwable normalFailure=invokeCloseOwnedInput(method, normal, null);
		check(normalFailure==null && normal.closes==1,
				"normal close state differed: closes="+normal.closes
				+" failure="+exceptionSummary(normalFailure));

		final IOException closeOnly=new IOException("synthetic close-only failure");
		final CloseTrackingInputStream closeOnlyStream=
				new CloseTrackingInputStream(closeOnly);
		final Throwable closeOnlyFailure=invokeCloseOwnedInput(method, closeOnlyStream, null);
		check(closeOnlyFailure==closeOnly && closeOnlyStream.closes==1,
				"close-only IOException was not thrown by identity: closes="
				+closeOnlyStream.closes+" failure="+exceptionSummary(closeOnlyFailure));

		final IOException sameIo=new IOException("synthetic same IO");
		final CloseTrackingInputStream sameIoStream=new CloseTrackingInputStream(sameIo);
		final Throwable sameIoFailure=invokeCloseOwnedInput(method, sameIoStream, sameIo);
		check(sameIoFailure==null && sameIoStream.closes==1
				&& sameIo.getSuppressed().length==0,
				"same IOException suppression differed: closes="+sameIoStream.closes
				+" failure="+exceptionSummary(sameIoFailure)
				+" suppressed="+exceptionSummaries(sameIo));

		final IOException primaryIo=new IOException("synthetic primary IO");
		final IOException closeAfterIo=new IOException("synthetic close after IO");
		final CloseTrackingInputStream ioStream=new CloseTrackingInputStream(closeAfterIo);
		final Throwable ioFailure=invokeCloseOwnedInput(method, ioStream, primaryIo);
		check(ioFailure==null && ioStream.closes==1
				&& hasOnlySuppressed(primaryIo, closeAfterIo),
				"IO primary suppression differed: closes="+ioStream.closes
				+" failure="+exceptionSummary(ioFailure)
				+" suppressed="+exceptionSummaries(primaryIo));

		final RuntimeException primaryRuntime=
				new RuntimeException("synthetic primary runtime");
		final IOException closeAfterRuntime=
				new IOException("synthetic close after runtime");
		final CloseTrackingInputStream runtimeStream=
				new CloseTrackingInputStream(closeAfterRuntime);
		final Throwable runtimeFailure=
				invokeCloseOwnedInput(method, runtimeStream, primaryRuntime);
		check(runtimeFailure==null && runtimeStream.closes==1
				&& hasOnlySuppressed(primaryRuntime, closeAfterRuntime),
				"RuntimeException primary suppression differed: closes="
				+runtimeStream.closes+" failure="+exceptionSummary(runtimeFailure)
				+" suppressed="+exceptionSummaries(primaryRuntime));

		final Error primaryError=new Error("synthetic primary error");
		final IOException closeAfterError=new IOException("synthetic close after error");
		final CloseTrackingInputStream errorStream=
				new CloseTrackingInputStream(closeAfterError);
		final Throwable errorFailure=invokeCloseOwnedInput(method, errorStream, primaryError);
		check(errorFailure==null && errorStream.closes==1
				&& hasOnlySuppressed(primaryError, closeAfterError),
				"Error primary suppression differed: closes="+errorStream.closes
				+" failure="+exceptionSummary(errorFailure)
				+" suppressed="+exceptionSummaries(primaryError));
	}

	private static Throwable invokeCloseOwnedInput(final Method method,
			final CloseTrackingInputStream input, final Throwable primary) throws Exception{
		try{
			method.invoke(null, input, primary);
			return null;
		}catch(InvocationTargetException e){
			return e.getCause();
		}
	}

	private static boolean hasOnlySuppressed(final Throwable primary,
			final Throwable suppressed){
		final Throwable[] array=primary.getSuppressed();
		return array.length==1 && array[0]==suppressed;
	}

	private static RuntimeException expectRuntime(final ThrowingCall call,
			final String message){
		try{
			call.run();
			throw new AssertionError(message);
		}catch(RuntimeException e){
			return e;
		}catch(Exception e){
			throw new AssertionError("Unexpected checked failure", e);
		}
	}

	private static boolean causeMessageContains(final RuntimeException failure,
			final String message){
		final Throwable cause=failure.getCause();
		return cause instanceof IOException && cause.getMessage()!=null
				&& cause.getMessage().contains(message);
	}

	private static boolean hasSuppressedMessage(final Throwable primary,
			final String message){
		for(final Throwable suppressed : primary.getSuppressed()){
			if(suppressed instanceof IOException && suppressed.getMessage()!=null
					&& suppressed.getMessage().contains(message)){
				return true;
			}
		}
		return false;
	}

	private static MainResult captureMain(final TrackingInputStream input,
			final String... args) throws Exception{
		final InputStream oldIn=System.in;
		final PrintStream oldOut=System.out;
		final PrintStream oldErr=System.err;
		final int oldFormat=CountGC.FORMAT;
		final boolean oldSummary=CountGC.SUMMARY_ONLY;
		final PrintStream oldOutstream=(PrintStream)readStaticField(CountGC.class, "outstream");
		final long oldLimsum=((Long)readStaticField(CountGC.class, "LIMSUM")).longValue();
		final int oldThreads=Shared.threads();
		final int oldBuffers=Shared.READ_BUFFER_NUM_BUFFERS();
		final boolean oldUseUnpigz=ReadWrite.USE_UNPIGZ;
		final boolean oldUseUnbgzip=ReadWrite.USE_UNBGZIP;
		final boolean oldPreferUnbgzip=ReadWrite.PREFER_UNBGZIP;
		final boolean oldUseGunzip=ReadWrite.USE_GUNZIP;
		final boolean oldAllowNativeBgzf=ReadWrite.ALLOW_NATIVE_BGZF;
		final boolean oldPreferNativeBgzfIn=ReadWrite.PREFER_NATIVE_BGZF_IN;
		final boolean oldPreferNativeBgzfOut=ReadWrite.PREFER_NATIVE_BGZF_OUT;
		final boolean oldRawmode=ReadWrite.RAWMODE;
		final ByteArrayOutputStream stdoutBytes=new ByteArrayOutputStream();
		final ByteArrayOutputStream stderrBytes=new ByteArrayOutputStream();
		final PrintStream stdout=new PrintStream(stdoutBytes, true,
				StandardCharsets.UTF_8.name());
		final PrintStream stderr=new PrintStream(stderrBytes, true,
				StandardCharsets.UTF_8.name());
		Throwable failure=null;
		MainSnapshot snapshot=null;
		try{
			System.setIn(input);
			System.setOut(stdout);
			System.setErr(stderr);
			ReadWrite.RAWMODE=false;
			CountGC.main(args);
		}catch(Throwable t){
			failure=t;
		} finally{
			snapshot=captureMainSnapshot();
			System.setIn(oldIn);
			System.setOut(oldOut);
			System.setErr(oldErr);
			CountGC.FORMAT=oldFormat;
			CountGC.SUMMARY_ONLY=oldSummary;
			setStaticField(CountGC.class, "outstream", oldOutstream);
			setStaticField(CountGC.class, "LIMSUM", Long.valueOf(oldLimsum));
			Shared.setThreads(oldThreads);
			Shared.setBuffers(oldBuffers);
			ReadWrite.USE_UNPIGZ=oldUseUnpigz;
			ReadWrite.USE_UNBGZIP=oldUseUnbgzip;
			ReadWrite.PREFER_UNBGZIP=oldPreferUnbgzip;
			ReadWrite.USE_GUNZIP=oldUseGunzip;
			ReadWrite.ALLOW_NATIVE_BGZF=oldAllowNativeBgzf;
			ReadWrite.PREFER_NATIVE_BGZF_IN=oldPreferNativeBgzfIn;
			ReadWrite.PREFER_NATIVE_BGZF_OUT=oldPreferNativeBgzfOut;
			ReadWrite.RAWMODE=oldRawmode;
			stdout.flush();
			stderr.flush();
			stdout.close();
			stderr.close();
		}
		return new MainResult(stdoutBytes.toString(StandardCharsets.UTF_8.name()),
				stderrBytes.toString(StandardCharsets.UTF_8.name()), failure, snapshot);
	}

	private static MainSnapshot captureMainSnapshot(){
		return new MainSnapshot(Shared.threads(), Shared.READ_BUFFER_NUM_BUFFERS(),
				ReadWrite.USE_UNPIGZ, ReadWrite.USE_UNBGZIP, ReadWrite.PREFER_UNBGZIP,
				ReadWrite.USE_GUNZIP, ReadWrite.ALLOW_NATIVE_BGZF,
				ReadWrite.PREFER_NATIVE_BGZF_IN, ReadWrite.PREFER_NATIVE_BGZF_OUT,
				ReadWrite.RAWMODE);
	}

	private static void requireDevFull(){
		final Path devFull=Paths.get(DEV_FULL);
		if(!Files.exists(devFull) || !Files.isWritable(devFull)){
			throw new AssertionError(DEV_FULL+" is required for outputerror scope");
		}
	}

	private static void checkOutputErrorDevice(final String label, final boolean baseline,
			final long[] expectedCounts, final WriterCountCall call) throws Exception{
		final WriterResult actual=captureDeviceOutputError(call);
		recordWriterObservation(label, actual, "outputPath", DEV_FULL);
		throwCleanupFailures(label, actual);
		check(actual.writers.isEmpty(),
				label+" output path left "+actual.writers.size()+" live writers");
		if(baseline){
			check(actual.failure==null && Arrays.equals(expectedCounts, actual.counts),
					label+" baseline failure/counts differed: "
					+exceptionSummary(actual.failure)+" "+Arrays.toString(actual.counts));
		}else{
			checkWriterOutputFailure(label, actual.failure);
		}
	}

	private static void checkOutputErrorPrimary(final String label, final boolean baseline,
			final String parserMessage, final IOException expectedReadFailure,
			final WriterCountCall call) throws Exception{
		final WriterResult actual=captureDeviceOutputError(call);
		recordWriterObservation(label, actual, "outputPath", DEV_FULL);
		throwCleanupFailures(label, actual);
		check(actual.writers.isEmpty(),
				label+" output path left "+actual.writers.size()+" live writers");
		if(expectedReadFailure!=null){
			check(actual.failure==expectedReadFailure,
					label+" IOException identity differed: "+exceptionSummary(actual.failure));
		}else{
			check(actual.failure instanceof RuntimeException && actual.failure.getMessage()!=null
					&& actual.failure.getMessage().contains(parserMessage),
					label+" parser diagnostic differed: "+exceptionSummary(actual.failure));
		}
		checkSuppressedOutputFailure(label, actual.failure, baseline ? 0 : 1);
	}

	private static WriterResult captureDeviceOutputError(final WriterCountCall call)
			throws Exception{
		final Set<TextStreamWriter> before=liveTextStreamWriters();
		final WriterResult actual=captureWriterCount(before, Paths.get(DEV_FULL), call);
		try{
			return actual;
		} finally{
			actual.cleanup();
		}
	}

	private static void checkWriterOutputFailure(final String label, final Throwable failure){
		check(isWriterOutputFailure(failure),
				label+" writer failure differed: "+exceptionSummary(failure));
	}

	private static void checkSuppressedOutputFailure(final String label,
			final Throwable primary, final int expected){
		final Throwable[] suppressed=(primary==null ? new Throwable[0] : primary.getSuppressed());
		check(suppressed.length==expected
				&& (expected==0 || isWriterOutputFailure(suppressed[0])),
				label+" suppressed writer failure differed: "
				+Arrays.toString(suppressed));
	}

	private static boolean isWriterOutputFailure(final Throwable failure){
		return failure instanceof IOException
				&& failure.getMessage()!=null
				&& failure.getMessage().contains("CountGC")
				&& failure.getMessage().contains(DEV_FULL);
	}

	private static void checkWriterCleanup(final Path scratch, final boolean baseline)
			throws Exception{
		final String oneRow="one\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";
		final String twoRow="two\t4\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.00000\t0.50000\n";
		final String okRow="ok\t5\t0.25000\t0.25000\t0.25000\t0.25000"
				+"\t0.20000\t0.50000\n";

		checkWriterNormal("normal FASTA writer", scratch.resolve("writer-normal.fa.tsv"),
				oneRow+twoRow, new long[]{2, 2, 2, 2, 1, 1}, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=">one\nACGTN\n>two\nACGT".getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFasta(new ByteArrayInputStream(input), out);
			}
		});

		checkWriterNormal("normal FASTQ writer", scratch.resolve("writer-normal.fq.tsv"),
				oneRow+twoRow, new long[]{2, 2, 2, 2, 1, 0}, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=("@one\nACGTN\n+\nIIIII\n"
						+"@two\nACGT\n+\nIIII\n").getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFastq(new ByteArrayInputStream(input), out);
			}
		});

		checkWriterFailure("parser FASTA writer", baseline,
				scratch.resolve("writer-parser.fa.tsv"), okRow,
				"Invalid FASTA sequence byte 0xC3", null, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=new byte[]{'>', 'o', 'k', '\n', 'A', 'C', 'G', 'T', 'N', '\n',
						'>', 'b', 'a', 'd', '\n', 'A', 'C', (byte)0xC3, 'T', '\n'};
				return CountGC.countFasta(new ByteArrayInputStream(input), out);
			}
		});

		checkWriterFailure("parser FASTQ writer", baseline,
				scratch.resolve("writer-parser.fq.tsv"), okRow,
				"Invalid FASTQ sequence byte 0xC3", null, new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=new byte[]{'@', 'o', 'k', '\n', 'A', 'C', 'G', 'T', 'N', '\n',
						'+', '\n', 'I', 'I', 'I', 'I', 'I', '\n', '@', 'b', 'a', 'd', '\n',
						'A', 'C', (byte)0xC3, 'T', '\n', '+', '\n', 'I', 'I', 'I', 'I', '\n'};
				return CountGC.countFastq(new ByteArrayInputStream(input), out);
			}
		});

		final IOException fastaFailure=new IOException("synthetic FASTA read failure");
		checkWriterFailure("read FASTA writer", baseline,
				scratch.resolve("writer-read.fa.tsv"), okRow, null, fastaFailure,
				new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=">ok\nACGTN\n>bad\nAC".getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFasta(new ThrowingInputStream(input, fastaFailure), out);
			}
		});

		final IOException fastqFailure=new IOException("synthetic FASTQ read failure");
		checkWriterFailure("read FASTQ writer", baseline,
				scratch.resolve("writer-read.fq.tsv"), okRow, null, fastqFailure,
				new WriterCountCall(){
			@Override
			public long[] count(final String out) throws Exception{
				final byte[] input=("@ok\nACGTN\n+\nIIIII\n"
						+"@bad\nAC").getBytes(StandardCharsets.US_ASCII);
				return CountGC.countFastq(new ThrowingInputStream(input, fastqFailure), out);
			}
		});
	}

	private static void checkWriterNormal(final String label, final Path output,
			final String expectedOutput, final long[] expectedCounts, final WriterCountCall call)
			throws Exception{
		final Set<TextStreamWriter> before=liveTextStreamWriters();
		Files.deleteIfExists(output);
		final WriterResult actual=captureWriterCount(before, output, call);
		final String actualOutput;
		try{
			actualOutput=readStringIfExists(output);
		} finally{
			actual.cleanup();
		}
		recordWriterObservation(label, actual, "output", actualOutput);
		throwCleanupFailures(label, actual);
		if(actual.failure instanceof Error){throw (Error)actual.failure;}
		if(actual.failure!=null){throw (Exception)actual.failure;}
		check(actual.writers.isEmpty(),
				label+" normal path left "+actual.writers.size()+" live writers");
		check(Arrays.equals(expectedCounts, actual.counts),
				label+" aggregate counts were "+Arrays.toString(actual.counts));
		check(expectedOutput.equals(actualOutput), label+" output differed");
	}

	private static void checkWriterFailure(final String label, final boolean baseline,
			final Path output, final String expectedOutput, final String parserMessage,
			final IOException expectedReadFailure, final WriterCountCall call) throws Exception{
		final Set<TextStreamWriter> before=liveTextStreamWriters();
		Files.deleteIfExists(output);
		final WriterResult actual=captureWriterCount(before, output, call);
		boolean started=false;
		boolean open=false;
		boolean hasPoison=true;
		String bufferText=null;
		try{
			if(actual.writers.size()==1){
				final TextStreamWriter writer=actual.writers.iterator().next();
				started=readBoolean(writer, "started");
				open=readBoolean(writer, "open");
				hasPoison=queueHasPoison(writer);
				bufferText=bufferText(writer);
			}
		} finally{
			actual.cleanup();
		}

		final String actualOutput=readStringIfExists(output);
		recordWriterObservation(label, actual, "output", actualOutput);
		throwCleanupFailures(label, actual);
		if(expectedReadFailure!=null){
			check(actual.failure==expectedReadFailure,
					label+" IOException identity differed: "+exceptionSummary(actual.failure));
		}else{
			check(actual.failure instanceof RuntimeException && actual.failure.getMessage()!=null
					&& actual.failure.getMessage().contains(parserMessage),
					label+" parser diagnostic differed: "+exceptionSummary(actual.failure));
		}
		if(baseline){
			check(actual.writers.size()==1, label+" new writer count was "
					+actual.writers.size());
			check(started, label+" writer was not marked started");
			check(open, label+" writer was not open before fixture cleanup");
			check(!hasPoison, label+" queue already held CountGC poison");
			check(expectedOutput.equals(bufferText),
					label+" retained buffer differed: "+bufferText);
			for(final Boolean cleanupState : actual.cleanupStates){
				check(!cleanupState.booleanValue(), label+" cleanup returned errorState=true");
			}
		}else{
			check(actual.writers.isEmpty(),
					label+" fixed path left "+actual.writers.size()+" live writers");
		}
		check(expectedOutput.equals(actualOutput), label+" post-cleanup output differed");
		check(newTextStreamWriters(before).isEmpty(),
				label+" left live writers after fixture cleanup");
	}

	private static Set<TextStreamWriter> newTextStreamWriters(final Set<TextStreamWriter> before){
		final Set<TextStreamWriter> current=liveTextStreamWriters();
		current.removeAll(before);
		return current;
	}

	private static Set<TextStreamWriter> liveTextStreamWriters(){
		final Set<TextStreamWriter> set=new HashSet<TextStreamWriter>();
		for(final Thread thread : Thread.getAllStackTraces().keySet()){
			if(thread instanceof TextStreamWriter){
				set.add((TextStreamWriter)thread);
			}
		}
		return set;
	}

	private static boolean readBoolean(final TextStreamWriter writer, final String field)
			throws Exception{
		return ((Boolean)readField(writer, field)).booleanValue();
	}

	private static String bufferText(final TextStreamWriter writer) throws Exception{
		@SuppressWarnings("unchecked")
		final ArrayList<CharSequence> buffer=
				(ArrayList<CharSequence>)readField(writer, "buffer");
		if(buffer==null){return null;}
		final StringBuilder sb=new StringBuilder();
		for(final CharSequence entry : buffer){sb.append(entry);}
		return sb.toString();
	}

	/** Captures new writer threads; the caller must inspect and clean them up. */
	private static WriterResult captureWriterCount(final Set<TextStreamWriter> before,
			final Path output, final WriterCountCall call) throws Exception{
		final int oldFormat=CountGC.FORMAT;
		long[] counts=null;
		Throwable failure=null;
		Set<TextStreamWriter> writers=new HashSet<TextStreamWriter>();
		try{
			CountGC.FORMAT=1;
			counts=call.count(output.toString());
		}catch(Throwable t){
			failure=t;
		} finally{
			CountGC.FORMAT=oldFormat;
			writers=newTextStreamWriters(before);
		}
		return new WriterResult(counts, failure, writers);
	}

	private static void recordWriterObservation(final String label,
			final WriterResult actual, final String outputLabel, final String output){
		System.out.println("OBS "+label
				+"\t"+outputLabel+"="+escape(output)
				+"\tcounts="+Arrays.toString(actual.counts)
				+"\texception="+exceptionSummary(actual.failure)
				+"\tsuppressed="+exceptionSummaries(actual.failure)
				+"\twriters="+actual.writers.size()
				+"\tcleanupStates="+actual.cleanupStates
				+"\tcleanupFailures="+actual.cleanupFailures.size());
	}

	private static void throwCleanupFailures(final String label,
			final WriterResult actual){
		if(actual.cleanupFailures.isEmpty()){return;}
		final AssertionError error=new AssertionError(label+" fixture cleanup failed");
		for(final Throwable failure : actual.cleanupFailures){
			error.addSuppressed(failure);
		}
		throw error;
	}

	private static boolean queueHasPoison(final TextStreamWriter writer) throws Exception{
		final Object poison=readStaticField(TextStreamWriter.class, "POISON2");
		for(final Object queued : (Iterable<?>)readField(writer, "queue")){
			if(queued==poison){return true;}
		}
		return false;
	}

	private static Object readField(final Object target, final String name) throws Exception{
		final Field field=target.getClass().getDeclaredField(name);
		field.setAccessible(true);
		return field.get(target);
	}

	private static Object readStaticField(final Class<?> clazz, final String name)
			throws Exception{
		final Field field=clazz.getDeclaredField(name);
		field.setAccessible(true);
		return field.get(null);
	}

	private static void setStaticField(final Class<?> clazz, final String name,
			final Object value) throws Exception{
		final Field field=clazz.getDeclaredField(name);
		field.setAccessible(true);
		field.set(null, value);
	}

	private static void checkArrayIndexFailure(final String label, final int expectedIndex,
			final CountCall call) throws Exception{
		final Captured actual=captureCount(call);
		recordObservation(label, actual);
		final StackTraceElement top=actual.failure==null || actual.failure.getStackTrace().length<1
				? null : actual.failure.getStackTrace()[0];
		final String message=actual.failure==null ? null : actual.failure.getMessage();
		check(actual.failure instanceof ArrayIndexOutOfBoundsException
				&& message!=null
				&& message.contains(Integer.toString(expectedIndex))
				&& top!=null
				&& "jgi.CountGC".equals(top.getClassName())
				&& "countFasta".equals(top.getMethodName()),
				label+" exception differed: "+exceptionSummary(actual.failure));
	}

	private static void checkClearInvalidSequence(final String label, final String message,
			final CountCall call) throws Exception{
		final Captured actual=captureCount(call);
		recordObservation(label, actual);
		if(actual.failure==null){
			throw new AssertionError(label+" unexpectedly completed");
		}
		final String actualMessage=actual.failure.getMessage();
		check(actual.failure instanceof RuntimeException
				&& actualMessage!=null && actualMessage.contains(message),
				label+" diagnostic differed: "+exceptionSummary(actual.failure));
	}

	private static void checkCapture(final String label, final long[] expectedCounts,
			final String expectedStdout, final CountCall call) throws Exception{
		final Captured actual=captureCount(call);
		recordObservation(label, actual);
		if(actual.failure!=null){throw actual.failure;}
		check(Arrays.equals(expectedCounts, actual.counts),
				label+" aggregate counts were "+Arrays.toString(actual.counts));
		check(expectedStdout.equals(actual.stdout), label+" stdout differed: "+actual.stdout);
	}

	private static Captured captureCount(final CountCall call) throws Exception{
		final int oldFormat=CountGC.FORMAT;
		final boolean oldSummary=CountGC.SUMMARY_ONLY;
		final PrintStream oldOut=System.out;
		final ByteArrayOutputStream capturedBytes=new ByteArrayOutputStream();
		final PrintStream captured=new PrintStream(capturedBytes, true,
				StandardCharsets.UTF_8.name());
		long[] counts=null;
		Exception failure=null;
		try{
			CountGC.FORMAT=1;
			CountGC.SUMMARY_ONLY=false;
			System.setOut(captured);
			counts=call.count();
			captured.flush();
		}catch(Exception e){
			captured.flush();
			failure=e;
		} finally{
			System.setOut(oldOut);
			CountGC.FORMAT=oldFormat;
			CountGC.SUMMARY_ONLY=oldSummary;
			captured.close();
		}
		return new Captured(capturedBytes.toString(StandardCharsets.UTF_8.name()),
				counts, failure);
	}

	private static void recordObservation(final String label, final Captured captured){
		System.out.println("OBS "+label
				+"\tstdout="+escape(captured.stdout)
				+"\tcounts="+Arrays.toString(captured.counts)
				+"\texception="+exceptionSummary(captured.failure));
	}

	private static String exceptionSummary(final Throwable exception){
		if(exception==null){return "none";}
		final StackTraceElement[] stack=exception.getStackTrace();
		final String top=stack.length>0 ? stack[0].toString() : "no-stack";
		return exception.getClass().getName()+": "+exception.getMessage()+" @ "+top;
	}

	private static String exceptionSummaries(final Throwable exception){
		final Throwable[] suppressed=(exception==null ? new Throwable[0] : exception.getSuppressed());
		final StringBuilder sb=new StringBuilder("[");
		for(int i=0; i<suppressed.length; i++){
			if(i>0){sb.append(", ");}
			sb.append(exceptionSummary(suppressed[i]));
		}
		return sb.append(']').toString();
	}

	private static String escape(final String text){
		return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n");
	}

	private static void checkOutEqualsPath(final Path scratch, final boolean baseline) throws Exception{
		final Path input=scratch.resolve("eq.fa");
		final Path outdir=Files.createDirectory(scratch.resolve("out=dir"));
		final Path target=outdir.resolve("eq.tsv");
		final Path truncated=scratch.resolve("out");
		final Path stdout=scratch.resolve("countgc.stdout");
		final Path stderr=scratch.resolve("countgc.stderr");
		Files.write(input, ">eqcase\nAACCGGTTNN\n".getBytes(StandardCharsets.US_ASCII));

		final Path countgc=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.exists(countgc), "countgc.sh not found at "+countgc);
		final int exit=run(scratch, "countgc", stdout, stderr, "bash", countgc.toString(),
				"-Xmx128m", "-Xms32m", "in="+input.toString(), "out="+target.toString(),
				"format=1", "ow=t", "simd=f");
		Files.write(scratch.resolve("countgc.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);
		check(exit==0, "countgc.sh exit="+exit
				+"\nstdout:\n"+readString(stdout)+"\nstderr:\n"+readString(stderr));

		final String expected="eqcase\t10\t0.25000\t0.25000\t0.25000\t0.25000\t0.20000\t0.50000\n";
		if(baseline){
			check(!Files.exists(target), "baseline unexpectedly wrote complete out= target");
			check(expected.equals(readString(truncated)),
					"baseline truncated sibling content differed");
		}else{
			check(!Files.exists(truncated), "fixed parser created truncated sibling");
			check(expected.equals(readString(target)), "fixed target content differed");
		}
	}

	private static void checkSamePathInputOutput(final Path scratch, final boolean baseline)
			throws Exception{
		final byte[] originalBytes=">samepath\nAACCGGTTNN\n".getBytes(StandardCharsets.US_ASCII);
		final Path original=scratch.resolve("samepath.original.fa");
		final Path same=scratch.resolve("samepath.fa");
		final Path stdout=scratch.resolve("samepath.stdout");
		final Path stderr=scratch.resolve("samepath.stderr");
		Files.write(original, originalBytes);
		Files.write(same, originalBytes);

		final Path countgc=Paths.get("countgc.sh").toAbsolutePath().normalize();
		check(Files.exists(countgc), "countgc.sh not found at "+countgc);
		final int exit=run(scratch, "samepath", stdout, stderr, "bash", countgc.toString(),
				"-Xmx128m", "-Xms32m", "in="+same.toString(), "out="+same.toString(),
				"format=1", "ow=t", "simd=f");
		Files.write(scratch.resolve("samepath.exit"), Arrays.asList(Integer.toString(exit)),
				StandardCharsets.US_ASCII);
		check("".equals(readString(stdout)),
				"same-path countgc.sh wrote unexpected stdout: "+readString(stdout));
		check(Arrays.equals(originalBytes, Files.readAllBytes(original)),
				"preserved same-path original changed");

		if(baseline){
			check(exit==0, "baseline same-path exit="+exit
					+"\nstdout:\n"+readString(stdout)+"\nstderr:\n"+readString(stderr));
			check(Arrays.equals(new byte[0], Files.readAllBytes(same)),
					"baseline same-path input was not truncated to zero bytes");
		}else{
			check(exit!=0, "fixed same-path input/output alias unexpectedly exited 0");
			check(readString(stderr).contains("was specified multiple times"),
					"fixed same-path stderr lacked duplicate-path diagnostic: "
					+readString(stderr));
			check(Arrays.equals(originalBytes, Files.readAllBytes(same)),
					"fixed same-path input bytes differed from preserved original");
		}
	}

	private static void checkIntGcNumeratorOverflow(final boolean baseline) throws Exception{
		final int oldFormat=CountGC.FORMAT;
		try{
			CountGC.FORMAT=2;
			final Method method=CountGC.class.getDeclaredMethod("toString2",
					StringBuilder.class, int[].class);
			method.setAccessible(true);
			final String actual=(String)method.invoke(null, new StringBuilder("gcOverflow"),
					new int[]{0, 1100000000, 1100000000, 0, 0, 0});
			final String expected=baseline ? "gcOverflow\t-0.95226\n" : "gcOverflow\t1.00000\n";
			check(expected.equals(actual), "FORMAT=2 overflow output differed: "+actual);
		}catch(InvocationTargetException e){
			throw new AssertionError(e.getCause());
		} finally{
			CountGC.FORMAT=oldFormat;
		}
	}

	private static int run(final Path scratch, final String stem, final Path stdout, final Path stderr,
			final String... command) throws Exception{
		Files.write(scratch.resolve(stem+".args"), Arrays.asList(command), StandardCharsets.US_ASCII);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.directory(Paths.get(".").toAbsolutePath().normalize().toFile());
		pb.redirectOutput(stdout.toFile());
		pb.redirectError(stderr.toFile());
		final Process process=pb.start();
		return process.waitFor();
	}

	private static int runWithPathFirst(final Path scratch, final String stem,
			final Path stdout, final Path stderr, final Path pathFirst, final Path pigzMarker,
			final String... command) throws Exception{
		Files.write(scratch.resolve(stem+".args"), Arrays.asList(command),
				StandardCharsets.US_ASCII);
		final ProcessBuilder pb=new ProcessBuilder(command);
		pb.directory(Paths.get(".").toAbsolutePath().normalize().toFile());
		pb.redirectOutput(stdout.toFile());
		pb.redirectError(stderr.toFile());
		final Map<String, String> env=pb.environment();
		final String path=env.get("PATH");
		env.put("PATH", pathFirst.toString()+File.pathSeparator+(path==null ? "" : path));
		env.put("FAKE_PIGZ_MARKER", pigzMarker.toString());
		final Process process=pb.start();
		Throwable primary=null;
		try{
			process.getOutputStream().close();
			if(!process.waitFor(40, TimeUnit.SECONDS)){
				throw new AssertionError("Owned countgc.sh supervisor exceeded 40 seconds; "
						+"inspect "+stderr);
			}
			return process.exitValue();
		}catch(Exception e){
			primary=e;
			throw e;
		}catch(Error e){
			primary=e;
			throw e;
		}finally{
			try{finishOwnedChild(process);}
			catch(RuntimeException e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
			catch(Error e){if(primary==null){throw e;}else if(primary!=e){primary.addSuppressed(e);}}
		}
	}

	private static void finishOwnedChild(final Process process){
		if(!process.isAlive()){return;}
		process.destroy();
		if(waitOwnedChild(process, 8)){return;}
		process.destroyForcibly();
		if(!waitOwnedChild(process, 5)){
			throw new AssertionError("Owned countgc.sh supervisor did not terminate "
					+"after forced cleanup.");
		}
	}

	private static boolean waitOwnedChild(final Process process, final int seconds){
		final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(seconds);
		boolean interrupted=false;
		try{
			while(process.isAlive()){
				final long remaining=deadline-System.nanoTime();
				if(remaining<=0){return false;}
				try{return process.waitFor(remaining, TimeUnit.NANOSECONDS);}
				catch(InterruptedException e){interrupted=true;}
			}
			return true;
		}finally{
			if(interrupted){Thread.currentThread().interrupt();}
		}
	}

	private static String readString(final Path path) throws IOException{
		return new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
	}

	private static String readStringIfExists(final Path path){
		try{
			return Files.exists(path) ? readString(path) : "<missing>";
		}catch(IOException ioe){
			return "<read failure: "+exceptionSummary(ioe)+">";
		}
	}

	private static final class Captured{
		Captured(final String stdout_, final long[] counts_, final Exception failure_){
			stdout=stdout_;
			counts=counts_;
			failure=failure_;
		}

		final String stdout;
		final long[] counts;
		final Exception failure;
	}

	private static final class MainResult{
		MainResult(final String stdout_, final String stderr_, final Throwable failure_,
				final MainSnapshot snapshot_){
			stdout=stdout_;
			stderr=stderr_;
			failure=failure_;
			snapshot=snapshot_;
		}

		final String stdout;
		final String stderr;
		final Throwable failure;
		final MainSnapshot snapshot;
	}

	private static final class MainSnapshot{
		MainSnapshot(final int threads_, final int buffers_, final boolean useUnpigz_,
				final boolean useUnbgzip_, final boolean preferUnbgzip_,
				final boolean useGunzip_, final boolean allowNativeBgzf_,
				final boolean preferNativeBgzfIn_, final boolean preferNativeBgzfOut_,
				final boolean rawmode_){
			threads=threads_;
			buffers=buffers_;
			useUnpigz=useUnpigz_;
			useUnbgzip=useUnbgzip_;
			preferUnbgzip=preferUnbgzip_;
			useGunzip=useGunzip_;
			allowNativeBgzf=allowNativeBgzf_;
			preferNativeBgzfIn=preferNativeBgzfIn_;
			preferNativeBgzfOut=preferNativeBgzfOut_;
			rawmode=rawmode_;
		}

		boolean isJavaGzipState(){
			return threads==1 && buffers>=2 && !useUnpigz && !useUnbgzip
					&& !preferUnbgzip && !useGunzip && !allowNativeBgzf
					&& !preferNativeBgzfIn && !preferNativeBgzfOut && !rawmode;
		}

		@Override
		public String toString(){
			return "threads="+threads+", buffers="+buffers+", useUnpigz="+useUnpigz
					+", useUnbgzip="+useUnbgzip+", preferUnbgzip="+preferUnbgzip
					+", useGunzip="+useGunzip+", allowNativeBgzf="+allowNativeBgzf
					+", preferNativeBgzfIn="+preferNativeBgzfIn
					+", preferNativeBgzfOut="+preferNativeBgzfOut
					+", rawmode="+rawmode;
		}

		final int threads;
		final int buffers;
		final boolean useUnpigz;
		final boolean useUnbgzip;
		final boolean preferUnbgzip;
		final boolean useGunzip;
		final boolean allowNativeBgzf;
		final boolean preferNativeBgzfIn;
		final boolean preferNativeBgzfOut;
		final boolean rawmode;
	}


	/** CountGC callback bound to one fixture output path. */
	private interface WriterCountCall{
		long[] count(String out) throws Exception;
	}

	private static final String DEV_FULL="/dev/full";

	/** CountGC writer run state; cleanup attempts every writer and stores all failures. */
	private static final class WriterResult{
		WriterResult(final long[] counts_, final Throwable failure_,
				final Set<TextStreamWriter> writers_){
			counts=counts_;
			failure=failure_;
			writers=writers_;
		}

		void cleanup(){
			for(final TextStreamWriter writer : writers){
				try{
					cleanupStates.add(Boolean.valueOf(writer.poisonAndWait()));
				}catch(Throwable t){
					cleanupFailures.add(t);
				}
			}
		}

		final long[] counts;
		final Throwable failure;
		final Set<TextStreamWriter> writers;
		final ArrayList<Boolean> cleanupStates=new ArrayList<Boolean>();
		final ArrayList<Throwable> cleanupFailures=new ArrayList<Throwable>();
	}

	/** InputStream that caps positive bulk reads while retaining InputStream semantics. */
	private static final class ChunkedInputStream extends InputStream{
		ChunkedInputStream(final byte[] data_, final int chunkSize_){
			if(chunkSize_<1){throw new IllegalArgumentException("chunkSize must be positive");}
			data=data_;
			chunkSize=chunkSize_;
		}

		@Override
		public int read(){
			if(pos>=data.length){return -1;}
			return data[pos++]&0xFF;
		}

		@Override
		public int read(final byte[] b){
			return read(b, 0, b.length);
		}

		@Override
		public int read(final byte[] b, final int off, final int len){
			if(b==null){throw new NullPointerException();}
			if(off<0 || len<0 || len>b.length-off){throw new IndexOutOfBoundsException();}
			if(len==0){return 0;}
			if(pos>=data.length){return -1;}
			final int toCopy=Math.min(Math.min(len, chunkSize), data.length-pos);
			System.arraycopy(data, pos, b, off, toCopy);
			pos+=toCopy;
			return toCopy;
		}

		private final byte[] data;
		private final int chunkSize;
		private int pos=0;
	}


	/** InputStream that emits fixed fixture bytes and then throws a fixed failure. */
	private static final class ThrowingInputStream extends InputStream{
		ThrowingInputStream(final byte[] data_, final IOException failure_){
			data=data_;
			failure=failure_;
		}

		@Override
		public int read() throws IOException{
			if(pos<data.length){return data[pos++]&0xFF;}
			throw failure;
		}

		@Override
		public int read(final byte[] b) throws IOException{
			return read(b, 0, b.length);
		}

		@Override
		public int read(final byte[] b, final int off, final int len) throws IOException{
			if(len<1){return 0;}
			if(pos>=data.length){throw failure;}
			final int toCopy=Math.min(len, data.length-pos);
			System.arraycopy(data, pos, b, off, toCopy);
			pos+=toCopy;
			return toCopy;
		}

		private final byte[] data;
		private final IOException failure;
		private int pos=0;
	}

	/** Borrowed System.in fixture that can prove CountGC.main does not close it. */
	private static final class TrackingInputStream extends InputStream{
		TrackingInputStream(final String text, final IOException failure_){
			data=text.getBytes(StandardCharsets.US_ASCII);
			failure=failure_;
		}

		TrackingInputStream(final byte[] data_, final IOException failure_){
			data=Arrays.copyOf(data_, data_.length);
			failure=failure_;
		}

		@Override
		public int read() throws IOException{
			if(pos<data.length){return data[pos++]&0xFF;}
			if(failure!=null){throw failure;}
			return -1;
		}

		@Override
		public int read(final byte[] b) throws IOException{
			return read(b, 0, b.length);
		}

		@Override
		public int read(final byte[] b, final int off, final int len)
				throws IOException{
			if(b==null){throw new NullPointerException();}
			if(off<0 || len<0 || len>b.length-off){throw new IndexOutOfBoundsException();}
			if(len==0){return 0;}
			if(pos>=data.length){
				if(failure!=null){throw failure;}
				return -1;
			}
			final int toCopy=Math.min(len, data.length-pos);
			System.arraycopy(data, pos, b, off, toCopy);
			pos+=toCopy;
			return toCopy;
		}

		@Override
		public void close(){
			closed=true;
		}

		private final byte[] data;
		private final IOException failure;
		private int pos=0;
		boolean closed=false;
	}

	private static final class CloseTrackingInputStream extends InputStream{
		CloseTrackingInputStream(final IOException closeFailure_){
			closeFailure=closeFailure_;
		}

		@Override
		public int read(){
			return -1;
		}

		@Override
		public void close() throws IOException{
			closes++;
			if(closeFailure!=null){throw closeFailure;}
		}

		private final IOException closeFailure;
		int closes=0;
	}

	private static final class BulkTrackingInputStream extends InputStream{
		BulkTrackingInputStream(final byte[] data_){
			data=Arrays.copyOf(data_, data_.length);
		}

		@Override
		public int read(){
			singleReads++;
			return -1;
		}

		@Override
		public int read(final byte[] b){
			arrayReads++;
			return copy(b, 0, b.length);
		}

		@Override
		public int read(final byte[] b, final int off, final int len){
			sliceReads++;
			return copy(b, off, len);
		}

		private int copy(final byte[] b, final int off, final int len){
			final int toCopy=Math.min(len, data.length);
			System.arraycopy(data, 0, b, off, toCopy);
			return toCopy;
		}

		private final byte[] data;
		int arrayReads=0;
		int sliceReads=0;
		int singleReads=0;
	}

	private static final class FinishedProcess extends Process{
		FinishedProcess(final int exit_){
			exit=exit_;
		}

		@Override
		public OutputStream getOutputStream(){
			return new ByteArrayOutputStream();
		}

		@Override
		public InputStream getInputStream(){
			return new ByteArrayInputStream(new byte[0]);
		}

		@Override
		public InputStream getErrorStream(){
			return new ByteArrayInputStream(new byte[0]);
		}

		@Override
		public int waitFor(){
			return exit;
		}

		@Override
		public int exitValue(){
			return exit;
		}

		@Override
		public void destroy(){}

		private final int exit;
	}

	private interface CountCall{
		long[] count() throws Exception;
	}

	private interface ThrowingCall{
		void run() throws Exception;
	}

	private static void check(final boolean condition, final String message){
		checks++;
		if(!condition){
			throw new AssertionError(message);
		}
	}

	private static String phase;
	private static String outdir;
	private static String scope="default";
	private static int checks;
}
