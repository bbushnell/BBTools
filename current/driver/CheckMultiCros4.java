package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import stream.MultiCros4;
import stream.Read;
import structures.ByteBuilder;

/**
 * Small ordinary-input round-trip check for MultiCros4 destination completion.
 * Uses public routing operations and literal FASTQ content checks without timing hooks.
 * Each invocation requires a fresh output directory and leaves its outputs for inspection.
 * @author Shinobu
 * @date October 2, 2026
 */
public final class CheckMultiCros4{

	/** Parses flags, runs one routing case and validates every emitted record.
	 * @param args outdir=path, reads=count, length=bases, shortlength=bases, streams=count,
	 * bounce=t/f, minreads=count, threaded=t/f */
	public static void main(String[] args){
		Thread.setDefaultUncaughtExceptionHandler(new StopOnError());
		if(!CheckMultiCros4.class.desiredAssertionStatus()){
			throw new IllegalStateException("This output check requires assertions enabled with -ea");
		}
		final CheckMultiCros4 check=new CheckMultiCros4(args);
		check.process();
		Shared.closeStream(check.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Uses the standard pre-parser and delegates common options to Parser.
	 * @param args Command-line flags */
	private CheckMultiCros4(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final int equals=arg.indexOf('=');
			final String key=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			final String value=(equals<0 ? null : arg.substring(equals+1));
			if(key.equals("outdir")){directory=new File(value);}
			else if(key.equals("reads")){reads=Parse.parseIntKMG(value);}
			else if(key.equals("length")){length=Parse.parseIntKMG(value);}
			else if(key.equals("shortlength")){shortLength=Parse.parseIntKMG(value);}
			else if(key.equals("streams")){streams=Parse.parseIntKMG(value);}
			else if(key.equals("bounce")){bounce=Parse.parseBoolean(value);}
			else if(key.equals("minreads")){minReads=Parse.parseIntKMG(value);}
			else if(key.equals("threaded")){threaded=Parse.parseBoolean(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		if(shortLength==-1){shortLength=length;}
		assert(directory!=null && reads>=3 && length>0 && shortLength>0 && streams>0 && minReads>=0) : "A fresh directory, at least three reads, positive lengths/stream count and nonnegative minimum are required";
		if(directory.exists() || !directory.mkdirs()){
			throw new IllegalArgumentException("Check output directory must be new: "+directory);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Routes fresh unpaired roots, closes normally, then checks names, order and content. */
	private void process(){
		final byte[] longBases=bases(length), longQuality=filled(length, (byte)40), longEncoded=filled(length, (byte)'I');
		final byte[][] bases={longBases, shortLength==length ? longBases : bases(shortLength)};
		final byte[][] quality={longQuality, shortLength==length ? longQuality : filled(shortLength, (byte)40)};
		final byte[][] encodedQuality={longEncoded, shortLength==length ? longEncoded : filled(shortLength, (byte)'I')};
		final MultiCros4 output=new MultiCros4(new File(directory, "sample_%.fq").getPath(),
			null, false, false, false, false, FileFormat.FASTQ, threaded, streams);
		output.readsPerBuffer=1;
		output.bytesPerBuffer=Integer.MAX_VALUE;
		output.minReadsToDump=minReads;
		if(threaded){output.start();}
		for(int i=0; i<reads; i++){
			final int type=(i%3==0 ? 0 : 1);
			final Read read=new Read(bases[type].clone(), quality[type].clone(), "r"+i, i);
			if(threaded){
				read.obj=target(i);
				final ArrayList<Read> list=new ArrayList<Read>(1);
				list.add(read);
				output.add(list);
			}else{output.add(read, target(i));}
		}
		output.close();
		assert(!output.errorState()) : "Normal routing reported an output error";
		final File[] files=directory.listFiles();
		final int countB=(bounce ? (int)((reads+1L)/3) : 0), countA=reads-countB;
		final int expectedFiles=(countA>=minReads ? 1 : 0)+(bounce && countB>=minReads ? 1 : 0);
		assert(files!=null && files.length==expectedFiles) : "Only destinations meeting minreads should be created";
		long observed=0, observedBases=0;
		long expectedResidual=0, expectedResidualBases=0;
		for(int destination=0; destination<(bounce ? 2 : 1); destination++){
			final String name=NAMES[destination];
			if((destination==0 ? countA : countB)<minReads){
				for(int id=0; id<reads; id++){
					if(name.equals(target(id))){
						expectedResidual++;
						expectedResidualBases+=(id%3==0 ? length : shortLength);
					}
				}
				continue;
			}
			final ByteFile input=ByteFile.makeByteFile(new File(directory, "sample_"+name+".fq").getPath(), false);
			final ByteBuilder header=new ByteBuilder();
			for(int id=0; id<reads; id++){
				if(!name.equals(target(id))){continue;}
				final int type=(id%3==0 ? 0 : 1);
				header.clear().append("@r").append(id);
				assert(Arrays.equals(input.nextLine(), header.toBytes())) : "Missing or reordered record "+id+" in "+name;
				assert(Arrays.equals(input.nextLine(), bases[type])) : "Sequence changed for record "+id;
				assert(Arrays.equals(input.nextLine(), PLUS)) : "FASTQ separator changed for record "+id;
				assert(Arrays.equals(input.nextLine(), encodedQuality[type])) : "Quality changed for record "+id;
				observed++;
				observedBases+=bases[type].length;
			}
			assert(input.nextLine()==null) : "Unexpected additional records in "+name;
			final boolean error=input.close();
			assert(!error) : "Output verification reader reported an error";
		}
		assert(output.dumpResidual(null)==expectedResidual) : "Below-minimum roots must remain available for residual handling";
		assert(output.residualBases==expectedResidualBases) : "Residual bases must follow the mixed-length fixture";
		assert(observed+expectedResidual==reads) : "Every submitted root must be emitted or retained as a below-minimum residual";
		final long longReads=(reads+2L)/3;
		assert(observedBases+expectedResidualBases==longReads*length+(reads-longReads)*shortLength) : "Mixed-length fixture totals must match the three-record pattern";
		outstream.println("PASS reads="+observed+" bases="+observedBases+" streams="+streams+" bounce="+bounce+" shortlength="+shortLength+
			" minreads="+minReads+" residual="+expectedResidual+" threaded="+threaded);
	}

	/** Alternates a short A,B,A pattern or retains all roots in A for the resident control.
	 * @param id Input root index
	 * @return Stable destination name */
	private String target(final int id){return NAMES[bounce && id%3==1 ? 1 : 0];}

	/** Generates the literal repeated ACGT sequence used by the independent output oracle.
	 * @param length Positive sequence length
	 * @return Fresh sequence bytes */
	private static byte[] bases(final int length){
		assert(length>0) : "The normal fixture requires nonempty sequences";
		final byte[] result=new byte[length];
		for(int i=0; i<length; i++){result[i]=BASES[i&3];}
		return result;
	}

	/** Creates one immutable fixture template for qualities or their expected ASCII encoding.
	 * @param length Positive record length
	 * @param value Repeated byte value
	 * @return Fresh template array; submitted Read arrays are cloned separately */
	private static byte[] filled(final int length, final byte value){
		assert(length>0) : "Quality templates must match nonempty fixture sequences";
		final byte[] result=new byte[length];
		Arrays.fill(result, value);
		return result;
	}

	/** Reports an unexpected failure and terminates the check, including any helper threads. */
	private static final class StopOnError implements Thread.UncaughtExceptionHandler{
		/** Delegates to BBTools' process-level diagnostic termination.
		 * @param thread Thread reporting the failure
		 * @param failure Uncaught failure */
		@Override
		public void uncaughtException(Thread thread, Throwable failure){KillSwitch.exceptionKill(failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Fresh output directory and small invocation parameters. */
	private File directory;
	private int reads=48, length=32768, shortLength=-1, streams=2, minReads=0;
	private boolean bounce=true, threaded=false;
	private final PrintStream outstream;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final byte[] BASES={'A','C','G','T'}, PLUS={'+'};
	private static final String[] NAMES={"A", "B"};
}
