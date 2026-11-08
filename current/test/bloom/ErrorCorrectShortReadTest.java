package test.bloom;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bloom.ErrorCorrect;
import bloom.KCountArray;
import bloom.KCountArray7MTA;
import bloom.KmerCount6;
import dna.Data;
import shared.Shared;
import stream.ConcurrentGenericReadInputStream;
import stream.FASTQ;
import stream.FastqReadInputStream;
import stream.Read;
import structures.ListNum;

/** Real FASTQ/CRIS cases for symmetric short-read skipping and continued mate processing.
 * The observer retains real parsed Reads and seeds short-read error metadata before
 * passing batches unchanged to production. No counting or correction algorithm is replaced.
 * Run in a fresh JVM with outdir pointing to a new directory; retains inputs and diagnostics.
 * @author Brian Bushnell, Sucrose */
public final class ErrorCorrectShortReadTest{

	/** Runs six input shapes through detection, trusted counting and untrusted counting controls. */
	public static void main(final String[] args) throws IOException{
		if(args.length!=1 || !args[0].startsWith("outdir=")){throw new IllegalArgumentException("Usage: outdir=<new directory>");}
		final Path dir=Paths.get(args[0].substring(7)).toAbsolutePath();
		require(!Files.exists(dir), "Result directory must be new.");
		Files.createDirectories(dir);
		final int oldThreads=Shared.threads(), oldAdvance=ErrorCorrect.MIN_ADVANCE, oldBurst=ErrorCorrect.MAX_ERROR_BURST;
		final boolean oldTest=FASTQ.TEST_INTERLEAVED, oldForce=FASTQ.FORCE_INTERLEAVED, oldDetect=FASTQ.DETECT_QUALITY;
		final boolean oldInfo=ErrorCorrect.OUTPUT_INFO, oldVerbose=ErrorCorrect.verbose;
		final byte oldOffset=FASTQ.ASCII_OFFSET;
		int failures=0;
		try{
			Shared.setThreads(1);
			FASTQ.TEST_INTERLEAVED=FASTQ.FORCE_INTERLEAVED=FASTQ.DETECT_QUALITY=false;
			FASTQ.ASCII_OFFSET=33;
			ErrorCorrect.MIN_ADVANCE=1;
			ErrorCorrect.MAX_ERROR_BURST=3;
			ErrorCorrect.OUTPUT_INFO=ErrorCorrect.verbose=false;
			KCountArray7MTA.setSeed(20261105);
			final String[][] shapes={{"CC", null}, {"CCCCCA", null}, {"CC", "CCCCCA"},
					{"CCCCCA", "CC"}, {"CC", "CC"}, {"CCCCCA", "CCCCCA"}};
			for(int mode=0; mode<3; mode++){
				for(int shape=0; shape<shapes.length; shape++){
					failures+=exercise(dir, mode, shape, shapes[shape][0], shapes[shape][1]);
				}
			}
		}finally{
			Shared.setThreads(oldThreads);
			FASTQ.TEST_INTERLEAVED=oldTest; FASTQ.FORCE_INTERLEAVED=oldForce; FASTQ.DETECT_QUALITY=oldDetect;
			FASTQ.ASCII_OFFSET=oldOffset;
			ErrorCorrect.MIN_ADVANCE=oldAdvance; ErrorCorrect.MAX_ERROR_BURST=oldBurst;
			ErrorCorrect.OUTPUT_INFO=oldInfo; ErrorCorrect.verbose=oldVerbose;
		}
		if(failures>0){throw new AssertionError("Short-read stream contract failed "+failures+" of 18 cases.");}
		System.out.println("ERROR_CORRECT_SHORT_READ_TEST_PASS cases=18");
	}

	/** Uses genuine input threads and closes them even when the captured baseline throws. */
	private static int exercise(final Path dir, final int mode, final int shape, final String first, final String second) throws IOException{
		final String name=(mode==0 ? "detect" : mode==1 ? "trusted" : "untrusted")+"-"+shape;
		final Path in1=dir.resolve(name+"-1.fq"), in2=(second==null ? null : dir.resolve(name+"-2.fq"));
		writeRead(in1, first, "pair/1");
		if(in2!=null){writeRead(in2, second, "pair/2");}
		final KCountArray trusted=KCountArray.makeNew(2000003, 4, 3), output=KCountArray.makeNew(2000003, 4, 3);
		trusted.increment(341, 5);//CCCCC
		trusted.increment(682, 5);//GGGGG
		for(int key=0; key<1024; key++){
			require(trusted.read(key)==(key==341 || key==682 ? 5 : 0) && output.read(key)==0,
					"Literal five-base tables must be collision-free before the streamed case.");
		}
		final FastqReadInputStream source1=new FastqReadInputStream(in1.toString(), false);
		final FastqReadInputStream source2=(in2==null ? null : new FastqReadInputStream(in2.toString(), false));
		final ObservedStream input=new ObservedStream(source1, source2, first, second);
		final Path diagnostic=dir.resolve(name+".log");
		final PrintStream oldOut=Data.sysout;
		String exception="none";
		try(PrintStream log=new PrintStream(Files.newOutputStream(diagnostic), true, "UTF-8")){
			Data.sysout=log;
			input.start();
			try{
				if(mode==0){ErrorCorrect.detect(input, trusted, 5, 2, -1, null);}
				else{KmerCount6.count(input, 5, false, output, mode==1 ? trusted : null, 2, 1, true);}
			}catch(AssertionError | NullPointerException caught){
				exception=caught.getClass().getSimpleName();
				caught.printStackTrace(log);
			}finally{input.close();}
		}finally{Data.sysout=oldOut;}
		require(!input.errorState(), "The actual FASTQ reader must close without an input error.");
		final int expectedReads=(second==null ? 1 : 2);
		require(input.observed.size()==expectedReads, "The real stream must expose the complete literal pair or singleton.");
		boolean pass=exception.equals("none");
		int longReads=0;
		for(Read read : input.observed){
			final boolean shortRead=read.length()<5;
			if(!shortRead){longReads++;}
			final String expected=(shortRead ? "CC" : mode==0 ? "CCCCCC" : mode==1 ? "CCCCCN" : "CCCCCA");
			final byte[] expectedQual=new byte[expected.length()];
			Arrays.fill(expectedQual, (byte)30);
			if(!shortRead && mode<2){expectedQual[5]=(byte)(mode==0 ? 4 : 0);}
			final int expectedErrors=(shortRead ? 7 : mode==1 ? 1 : 0);
			final boolean readPass=Arrays.equals(read.bases, expected.getBytes(StandardCharsets.US_ASCII)) &&
					Arrays.equals(read.quality, expectedQual) && read.errors==expectedErrors;
			pass&=readPass;
			System.out.println(name+" read="+read.id+" expected="+expected+" actual="+new String(read.bases, StandardCharsets.US_ASCII)+
					" expectedErrors="+expectedErrors+" actualErrors="+read.errors+" pass="+readPass);
		}
		if(mode==0 && exception.equals("none")){
			final Matcher corrected=Pattern.compile("Reads Corrected:\\s+(\\d+)").matcher(new String(Files.readAllBytes(diagnostic), StandardCharsets.UTF_8));
			pass&=corrected.find() && Integer.parseInt(corrected.group(1))==longReads;
		}else if(mode!=0){
			pass&=output.read(341)==longReads && output.read(340)==(mode==2 ? longReads : 0);
		}
		System.out.println(name+" exception="+exception+" shortMetadataExpected=7 longReads="+longReads+" pass="+pass);
		return pass ? 0 : 1;
	}

	/** Writes one unambiguous ASCII-33 read with quality thirty. */
	private static void writeRead(final Path path, final String bases, final String id) throws IOException{
		final char[] qualities=new char[bases.length()];
		Arrays.fill(qualities, '?');
		Files.write(path, ("@"+id+"\n"+bases+"\n+\n"+new String(qualities)+"\n").getBytes(StandardCharsets.US_ASCII));
	}

	/** Distinguishes invalid fixture setup from the baseline's targeted AssertionError/NPE. */
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalStateException(message);}
	}

	/** Observes actual parsed records; delegates all lifecycle, batching and pairing to production CRIS. */
	private static final class ObservedStream extends ConcurrentGenericReadInputStream{
		ObservedStream(final FastqReadInputStream first, final FastqReadInputStream second, final String firstBases, final String secondBases){
			super(first, second, -1);
			this.firstBases=firstBases; this.secondBases=secondBases;
		}
		@Override
		public synchronized ListNum<Read> nextList(){
			final ListNum<Read> batch=super.nextList();
			if(batch!=null && batch.list!=null){
				for(Read root : batch.list){
					require(observed.isEmpty(), "Fixture contains exactly one input root.");
					observe(root, firstBases);
					require((root.mate!=null)==(secondBases!=null), "Actual stream pairing must match the literal fixture.");
					if(root.mate!=null){observe(root.mate, secondBases);}
				}
			}
			return batch;
		}
		private void observe(final Read read, final String expected){
			require(Arrays.equals(read.bases, expected.getBytes(StandardCharsets.US_ASCII)), "Actual parser must preserve the input sequence.");
			require(read.quality!=null && read.quality.length==expected.length(), "Actual parser must produce matching qualities.");
			for(byte quality : read.quality){require(quality==30, "Actual parser must decode ASCII-33 quality thirty.");}
			read.errors=(read.length()<5 ? 7 : 0);
			observed.add(read);
		}
		final ArrayList<Read> observed=new ArrayList<Read>();
		final String firstBases, secondBases;
	}
}
