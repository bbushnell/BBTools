package driver;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.KillSwitch;
import shared.Shared;
import stream.ArrayListSet;
import stream.ConcurrentReadOutputStream;
import stream.MultiCros;
import stream.Read;
import structures.ByteBuilder;
import template.ThreadWaiter;

/**
 * Checks ordinary concurrent submissions to an unordered MultiCros.
 * Producers own their lists, buckets and Read payloads; only the output is shared.
 * Checks registration identity and literal FASTQ records after producers and writers join.
 * Output order between producers is deliberately unspecified. This is a small
 * correctness check, not a race reproduction, scheduling test or benchmark.
 * @author Shinobu
 * @date October 2, 2026
 */
public final class CheckMultiCros{

	/** Runs one fresh-directory check with assertions enabled.
	 * @param args outdir=path, producers=1..4, destinations=1..4, grouped=t/f, paired=t/f */
	public static void main(String[] args){
		Thread.setDefaultUncaughtExceptionHandler(new StopOnError());
		if(!CheckMultiCros.class.desiredAssertionStatus()){
			throw new IllegalStateException("CheckMultiCros requires -ea for its output oracle");
		}
		final CheckMultiCros check=new CheckMultiCros(args);
		check.process();
		Shared.closeStream(check.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Expands configs and delegates shared flags to the standard Parser.
	 * @param args Command-line flags */
	private CheckMultiCros(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		outstream=pp.outstream;
		final Parser parser=new Parser();
		for(String arg : pp.args){
			final int equals=arg.indexOf('=');
			final String key=(equals<0 ? arg : arg.substring(0, equals)).toLowerCase(Locale.ROOT);
			final String value=(equals<0 ? null : arg.substring(equals+1));
			if(key.equals("outdir")){directory=new File(value);}
			else if(key.equals("producers")){producers=Parse.parseIntKMG(value);}
			else if(key.equals("destinations")){destinations=Parse.parseIntKMG(value);}
			else if(key.equals("grouped")){grouped=Parse.parseBoolean(value);}
			else if(key.equals("paired")){paired=Parse.parseBoolean(value);}
			else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown flag: "+arg);}
		}
		Parser.processQuality();
		assert(directory!=null && producers>=1 && producers<=4 && destinations>=1 && destinations<=4)
			: "Bounded normal fixture requires a directory, 1..4 producers and 1..4 destinations";
		if(!directory.mkdirs()){throw new IllegalArgumentException("Output directory must be fresh: "+directory);}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Runs producers, then closes children after submissions stop and verifies all output. */
	private void process(){
		final MultiCros output=new MultiCros(new File(directory, "sample_%.fq").getPath(),
			null, false, false, false, false, false, FileFormat.FASTQ, 8);
		final String[] names=new String[destinations];
		for(int d=0; d<names.length; d++){names[d]="D"+d;}
		final ArrayList<Producer> workers=new ArrayList<Producer>(producers);
		for(int p=0; p<producers; p++){workers.add(new Producer(output, names, p, paired, grouped));}
		final boolean joined=ThreadWaiter.startAndWait(workers);
		assert(joined) : "All producers must join before traversing MultiCros collections";
		for(Producer worker : workers){assert(worker.success) : "Producer did not finish its normal submissions";}
		final boolean error=ReadWrite.closeStreams(output);
		assert(!error && !output.errorState() && output.finishedSuccessfully()) : "Normal output must finish before verification";
		assert(output.streamMap.size()==destinations && output.streamList.size()==destinations)
			: "Each named destination must register exactly one child";
		for(int d=0; d<destinations; d++){
			final ConcurrentReadOutputStream child=output.streamMap.get(names[d]);
			assert(child!=null && output.streamList.contains(child)) : "Registered map/list children must agree";
			for(Producer worker : workers){
				assert(worker.children[d]==child) : "Every producer must observe the same child for "+names[d];
			}
		}
		verify(names);
		final int total=producers*ROOTS_PER_PRODUCER*(paired ? 2 : 1);
		outstream.println("PASS reads="+total+" bases="+(total*BASES.length)+" producers="+producers+
			" destinations="+destinations+" grouped="+grouped+" paired="+paired);
	}

	/** Verifies membership, uniqueness, destination and literal content without imposing order.
	 * @param names Exactly the destination names submitted by every producer */
	private void verify(final String[] names){
		final int mates=paired ? 2 : 1, roots=producers*ROOTS_PER_PRODUCER;
		final boolean[] seen=new boolean[roots*mates];
		final ByteBuilder expectedHeader=new ByteBuilder();
		final File[] files=directory.listFiles();
		assert(files!=null && files.length==destinations) : "Only one FASTQ file per destination should exist";
		int count=0;
		for(int d=0; d<names.length; d++){
			final ByteFile input=ByteFile.makeByteFile(new File(directory, "sample_"+names[d]+".fq").getPath(), false);
			for(byte[] header=input.nextLine(); header!=null; header=input.nextLine()){
				assert(header.length>=5 && header[0]=='@' && header[1]=='r' && header[header.length-2]=='/')
					: "Expected literal @r<id>/<mate> fixture header";
				final int id=Parse.parseInt(header, 2, header.length-2), mate=header[header.length-1]-'1';
				assert(id>=0 && id<roots && mate>=0 && mate<mates) : "Unexpected fixture read identifier or mate";
				expectedHeader.clear().append("@r").append(id).append('/').append(mate+1);
				assert(Arrays.equals(header, expectedHeader.toBytes())) : "Literal fixture header changed";
				assert((id%ROOTS_PER_PRODUCER/BATCH_SIZE)%destinations==d) : "Read routed to the wrong destination";
				final int index=id*mates+mate;
				assert(!seen[index]) : "Duplicate output for root="+id+" mate="+mate;
				seen[index]=true;
				assert(Arrays.equals(input.nextLine(), BASES)) : "Literal sequence changed for root="+id;
				assert(Arrays.equals(input.nextLine(), PLUS)) : "Literal FASTQ separator changed for root="+id;
				assert(Arrays.equals(input.nextLine(), ENCODED_QUALITY)) : "Literal phred33 quality changed for root="+id;
				count++;
			}
			final boolean error=input.close();
			assert(!error) : "Normal output verification reader reported an error";
		}
		assert(count==seen.length) : "Every submitted root and mate must appear exactly once: "+count+" != "+seen.length;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Owns all mutable input state for one ordinary producer. */
	private static final class Producer extends Thread{
		/** Captures shared output and immutable names; allocates private identity observations. */
		Producer(MultiCros output_, String[] names_, int number_, boolean paired_, boolean grouped_){
			output=output_; names=names_; number=number_; paired=paired_; grouped=grouped_;
			children=new ConcurrentReadOutputStream[names.length];
		}

		/** Submits independent batches, with per-thread buckets for the Seal-style path. */
		@Override
		public void run(){
			final ArrayListSet buckets=(grouped ? new ArrayListSet(false) : null);
			for(int start=0; start<ROOTS_PER_PRODUCER; start+=BATCH_SIZE){
				final int destination=(start/BATCH_SIZE)%names.length;
				final String name=names[destination];
				final ConcurrentReadOutputStream child=output.getStream(name);
				assert(children[destination]==null || children[destination]==child)
					: "Cached child identity must stay stable throughout ordinary submissions";
				children[destination]=child;
				final ArrayList<Read> list=(grouped ? null : new ArrayList<Read>(BATCH_SIZE));
				for(int i=start; i<start+BATCH_SIZE; i++){
					final int id=number*ROOTS_PER_PRODUCER+i;
					final Read first=read(id, 0);
					if(paired){
						final Read second=read(id, 1);
						first.mate=second; second.mate=first;
					}
					if(grouped){buckets.add(first, name);}
					else{list.add(first);}
				}
				final long listID=(long)number*(ROOTS_PER_PRODUCER/BATCH_SIZE)+start/BATCH_SIZE;
				if(grouped){output.add(buckets, listID);}
				else{output.add(list, listID, name);}
			}
			success=true;
		}

		/** Makes a fresh payload that remains borrowed by the writer until completion. */
		private static Read read(final int id, final int mate){
			assert(id>=0 && (mate==0 || mate==1)) : "Fixture identifiers and pair flags require a nonnegative root and mate0/1";
			final byte[] quality=new byte[BASES.length];
			Arrays.fill(quality, (byte)40);
			final Read read=new Read(BASES.clone(), quality, "r"+id+"/"+(mate+1), id);
			read.setPairnum(mate);
			return read;
		}

		final MultiCros output;
		final String[] names;
		final int number;
		final boolean paired, grouped;
		final ConcurrentReadOutputStream[] children;
		boolean success;
	}

	/** Ends this check if an unexpected producer/output error otherwise strands helper threads. */
	private static final class StopOnError implements Thread.UncaughtExceptionHandler{
		/** Reports the unexpected failure through BBTools' process-level termination helper. */
		@Override
		public void uncaughtException(Thread thread, Throwable failure){KillSwitch.exceptionKill(failure);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private File directory;
	private int producers=4, destinations=4;
	private boolean grouped=false, paired=false;
	private final PrintStream outstream;

	/*--------------------------------------------------------------*/
	/*----------------           Constants          ----------------*/
	/*--------------------------------------------------------------*/

	private static final int ROOTS_PER_PRODUCER=128, BATCH_SIZE=8;
	private static final byte[] BASES={'A','C','G','T','A','C','G','T'};
	private static final byte[] ENCODED_QUALITY={'I','I','I','I','I','I','I','I'}, PLUS={'+'};
}
