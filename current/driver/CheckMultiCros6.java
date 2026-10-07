package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import stream.BufferedMultiCros;
import stream.MultiCros6;
import stream.Read;
import structures.ByteBuilder;

/**
 * Checks ordinary MultiCros6 routing and its diagnostic main with a small valid fixture.
 * Direct and batch modes observe reported destinations before close; all modes verify
 * exact FASTQ contents after close. This is not a heap-footprint or performance measurement.
 * @author Shinobu
 * @date October 2, 2026
 */
public final class CheckMultiCros6{

	/** Runs one fresh-directory routing check.
	 * @param args outdir=path, mode=direct/batch/main, reads=count, length=bases,
	 * expectbefore=count, verbose=t/f, and standard Parser/BufferedMultiCros settings */
	public static void main(String[] args){
		Thread.setDefaultUncaughtExceptionHandler(new StopOnError());
		if(!CheckMultiCros6.class.desiredAssertionStatus()){
			throw new IllegalStateException("This check requires assertions enabled with -ea");
		}
		final CheckMultiCros6 check=new CheckMultiCros6(args);
		check.process();
		Shared.closeStream(check.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses fixture options and forwards ordinary settings to BBTools parsers.
	 * @param args Command-line flags */
	private CheckMultiCros6(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg : args){
			final int equals=arg.indexOf('=');
			final String key=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase();
			final String value=(equals<0 ? null : arg.substring(equals+1));
			if(key.equals("outdir")){directory=new File(value);}
			else if(key.equals("mode")){mode=value;}
			else if(key.equals("reads")){reads=Parse.parseIntKMG(value);}
			else if(key.equals("length")){length=Parse.parseIntKMG(value);}
			else if(key.equals("expectbefore")){expectedBefore=Parse.parseIntKMG(value);}
			else if(key.equals("verbose")){BufferedMultiCros.verbose=Parse.parseBoolean(value);}
			else if(!BufferedMultiCros.parseStatic(arg, key, value) && !parser.parse(arg, key, value)){
				throw new IllegalArgumentException("Unknown flag: "+arg);
			}
		}
		Parser.processQuality();
		assert(directory!=null && reads>=2 && length>0) : "The two-destination fixture needs a fresh directory and positive records";
		assert(mode.equals("direct") || mode.equals("batch") || mode.equals("main")) : "Unknown routing mode: "+mode;
		assert(mode.equals("main") || expectedBefore>=0) : "Specify the expected pre-close destination count for this routing check";
		if(directory.exists() || !directory.mkdirs()){
			throw new IllegalArgumentException("Check directory must be new: "+directory);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Creates fresh payloads, performs the selected routing path and checks every output byte. */
	private void process(){
		final byte[] bases=new byte[length], quality=new byte[length], encoded=new byte[length];
		for(int i=0; i<length; i++){bases[i]=BASES[i&3];}
		Arrays.fill(quality, (byte)40);
		Arrays.fill(encoded, (byte)'I');
		final File outputDirectory=new File(directory, "output");
		final boolean created=outputDirectory.mkdir();
		assert(created) : "Fresh output directory is needed to reject stale fixture data";
		final String pattern=new File(outputDirectory, "sample_%.fq").getPath();
		int before=-1;
		if(mode.equals("main")){
			final File input=new File(directory, "input.fq");
			final ByteStreamWriter writer=new ByteStreamWriter(input.getPath(), false, false, false);
			writer.start();
			final ByteBuilder bb=new ByteBuilder();
			for(int i=0; i<reads; i++){
				bb.clear().append('@').append(header(i)).nl().append(bases).nl().append('+').nl().append(encoded).nl();
				writer.print(bb);
			}
			final boolean error=writer.poisonAndWait();
			assert(!error) : "The valid diagnostic-main input must be complete before reading it";
			MultiCros6.main(new String[]{input.getPath(), pattern});
		}else{
			final MultiCros6 output=new MultiCros6(pattern, null, false, false, false, false, FileFormat.FASTQ, false, 4);
			output.minReadsToDump=0;
			for(int i=0; i<reads; i++){
				final Read read=new Read(bases.clone(), quality.clone(), header(i), i);
				if(mode.equals("direct")){output.add(read, NAMES[i&1]);}
				else{
					read.obj=NAMES[i&1];
					final ArrayList<Read> list=new ArrayList<Read>(1);
					list.add(read);
					output.add(list);
				}
			}
			before=0;
			for(byte b : output.report().toBytes()){if(b=='\n'){before++;}}
			assert(before==expectedBefore) : "Reported destinations before close: "+before+", expected "+expectedBefore;
			output.close();
			assert(!output.errorState()) : "Normal output close reported an error";
		}
		final File[] files=outputDirectory.listFiles();
		assert(files!=null && files.length==2) : "Exactly two barcode destinations must be created";
		int observed=0;
		final ByteBuilder bb=new ByteBuilder();
		for(int destination=0; destination<2; destination++){
			final ByteFile input=ByteFile.makeByteFile(new File(outputDirectory, "sample_"+NAMES[destination]+".fq").getPath(), false);
			for(long index=destination; index<reads; index+=2){
				final int id=(int)index;
				bb.clear().append('@').append(header(id));
				assert(Arrays.equals(input.nextLine(), bb.toBytes())) : "Missing or reordered root "+id;
				assert(Arrays.equals(input.nextLine(), bases)) : "Changed sequence for root "+id;
				assert(Arrays.equals(input.nextLine(), PLUS)) : "Changed separator for root "+id;
				assert(Arrays.equals(input.nextLine(), encoded)) : "Changed quality for root "+id;
				observed++;
			}
			assert(input.nextLine()==null) : "Unexpected extra records in "+NAMES[destination];
			final boolean error=input.close();
			assert(!error) : "Verification reader reported an error";
		}
		assert(observed==reads) : "Each submitted root must appear once after close";
		outstream.println("PASS mode="+mode+" before="+before+" reads="+observed+" bases="+((long)observed*length));
	}

	/** Builds a modern Illumina header with one of two literal barcode strings.
	 * @param id Unique root index
	 * @return Full header without its FASTQ at-sign */
	private static String header(final int id){return "fixture:1:flow:1:1:1:"+id+" 1:N:0:"+NAMES[id&1];}

	/** Terminates a failed check rather than leaving output helper threads alive. */
	private static final class StopOnError implements Thread.UncaughtExceptionHandler{
		/** Reports an unexpected failure through BBTools' standard process termination.
		 * @param thread Reporting thread
		 * @param failure Uncaught failure */
		@Override
		public void uncaughtException(Thread thread, Throwable failure){KillSwitch.exceptionKill(failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Fresh fixture root, routing selector and deterministic fixture dimensions. */
	private File directory;
	private String mode="direct";
	private int reads=12, length=32768, expectedBefore=-1;
	private final PrintStream outstream;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final String[] NAMES={"ACGT", "TGCA"};
	private static final byte[] BASES={'A','C','G','T'}, PLUS={'+'};
}
