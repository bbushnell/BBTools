package assemble;

import java.io.File;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Retains all training candidates and an exact, seeded uniform control subset.
 * Consumes a validated corpus and its row sidecar; never opens held-out tables.
 * @author Fischl
 */
public final class SameKTrainingSample {

	/** Parses standard config/JVM options and creates a new training-only directory. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SameKTrainingSample.class, false);
		try{new SameKTrainingSample(pp.args, pp.outstream).run();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Requires the master training table, full row sidecar and an explicit quota/seed. */
	private SameKTrainingSample(final String[] args, final PrintStream stream){
		outstream=stream;
		final Parser parser=new Parser();
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String name=split[0].toLowerCase(Locale.ROOT);
			while(name.startsWith("-")){name=name.substring(1);}
			final String value=split.length>1 ? split[1] : null;
			if(name.equals("rows")){rows=value;}
			else if(name.equals("controls")){quota=Integer.parseInt(value);}
			else if(name.equals("seed")){seed=Long.parseLong(value);}
			else if(!parser.parse(arg, name, value)){throw new IllegalArgumentException("Unknown sampling option: "+arg);}
		}
		input=parser.in1; output=parser.out1;
		if(input==null || rows==null || output==null || quota<0 || seed<0 || parser.in2!=null || parser.out2!=null ||
				!Tools.testInputFiles(false, true, input, rows) || new File(output).exists()){
			throw new IllegalArgumentException("Required in=train.tsv rows=rows.tsv out=NEW_directory controls=N seed=N.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts only training controls before streaming a without-replacement sample. */
	private void run(){
		final long[] totals=scan(false, null, null, 0);
		if(totals[1]>Integer.MAX_VALUE || quota>totals[1]){
			throw new IllegalArgumentException("Control quota exceeds population or supported integer range: "+quota+" / "+totals[1]);
		}
		if(!new File(output).mkdir()){throw new IllegalArgumentException("Output directory must be new: "+output);}
		final ByteStreamWriter table=open(output+"/train.tsv"), provenance=open(output+"/rows.tsv");
		long[] selected;
		boolean error=false;
		try{selected=scan(true, table, provenance, (int)totals[1]);}
		finally{error=table.poisonAndWait() | provenance.poisonAndWait();}
		if(error || !Arrays.equals(totals, selected)){
			throw new IllegalStateException("Input changed between passes or output failed; no completed subset exists.");
		}
		final String summary="SAME_K_TRAINING_SAMPLE_PASS candidates="+totals[0]+" controls="+quota+
				" available_controls="+totals[1]+" emitted="+(totals[0]+quota)+" seed="+seed;
		final ByteStreamWriter complete=open(output+"/complete.txt");
		complete.print(summary+"\n");
		if(complete.poisonAndWait()){throw new IllegalStateException("Cannot save subset completion marker.");}
		outstream.println(summary);
	}

	/** Validates row order/targets in both passes; feature bytes and source ordinals are preserved. */
	private long[] scan(final boolean emit, final ByteStreamWriter out, final ByteStreamWriter provenance, final int controls){
		assert(!emit || (out!=null && provenance!=null && controls>=quota)) : "Sampling requires an available quota and paired outputs.";
		final ByteFile meta=ByteFile.makeByteFile(rows, true), data=ByteFile.makeByteFile(input, true);
		final LineParser1 mp=new LineParser1('\t'), dp=new LineParser1('\t');
		final ByteBuilder buffer=new ByteBuilder();
		final long[] counts=new long[2], ordinals=new long[3];
		final Random random=new Random(seed);
		int remaining=controls, needed=quota;
		long emitted=0;
		boolean error=false;
		try{
			require(meta.nextLine(), "fold\trow\tgenome\tkind\tevent\ttarget");
			require(data.nextLine(), "#dims\t27\t1");
			final byte[] comment=data.nextLine();
			require(comment, "#schema=samek27-v1; unknown labels excluded; no class balancing");
			if(emit){
				out.print("#dims\t27\t1\n#schema=samek27-v1; training controls sampled uniformly without replacement\n");
				provenance.print("row\tsource_row\tgenome\tkind\tevent\ttarget\n");
			}
			for(byte[] line=meta.nextLine(); line!=null; line=meta.nextLine()){
				mp.set(line);
				if(mp.terms()!=6){throw new IllegalArgumentException("Wrong provenance width.");}
				final int fold=mp.termEquals("train", 0) ? 0 : mp.termEquals("validation", 0) ? 1 : mp.termEquals("test", 0) ? 2 : -1;
				final int kind=mp.termEquals("candidate", 3) ? 0 : mp.termEquals("control", 3) ? 1 : -1;
				final int target=mp.parseInt(5);
				if(fold<0 || kind<0 || mp.parseLong(1)!=++ordinals[fold] || mp.parseLong(4)<1 || (target!=0 && target!=1)){
					throw new IllegalArgumentException("Invalid fold/kind/ordinal/event/target in row sidecar.");
				}
				if(fold!=0){continue;}
				final byte[] vector=data.nextLine();
				if(vector==null){throw new IllegalArgumentException("Training table truncated before sidecar.");}
				dp.set(vector);
				if(dp.terms()!=28 || dp.parseInt(27)!=target){throw new IllegalArgumentException("Training width/target disagrees with sidecar.");}
				for(int i=0; i<27; i++){
					final float value=dp.parseFloat(i);
					if(!Float.isFinite(value) || value<0){throw new IllegalArgumentException("Invalid training feature.");}
				}
				counts[kind]++;
				if(!emit){continue;}
				boolean keep=true;
				if(kind==1){
					if(remaining<1){throw new IllegalArgumentException("Control population grew between passes.");}
					// Conditional inclusion probability selects exactly quota of the original controls.
					keep=random.nextInt(remaining)<needed;
					remaining--; if(keep){needed--;}
				}
				if(!keep){continue;}
				out.print(buffer.clear().append(vector).nl());
				int start=0;
				while(line[start]!='\t'){start++;}
				provenance.print(buffer.clear().append(++emitted).append(line, start, line.length-start).nl());
			}
			if(data.nextLine()!=null){throw new IllegalArgumentException("Training table has unpaired extra rows.");}
			if(emit && (remaining!=0 || needed!=0 || emitted!=counts[0]+quota)){
				throw new IllegalStateException("Subset count conservation failed.");
			}
		}finally{error=meta.close() | data.close();}
		if(error){throw new IllegalStateException("Incomplete sampling input read.");}
		return counts;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Static Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Starts a new output writer; callers observe its terminal error status. */
	private static ByteStreamWriter open(final String path){
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		return writer;
	}

	/** Exact headers prevent routing the wrong dataset format through a plausible parser. */
	private static void require(final byte[] actual, final String expected){
		if(!Arrays.equals(actual, expected.getBytes(StandardCharsets.UTF_8))){throw new IllegalArgumentException("Expected header: "+expected);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final PrintStream outstream;
	private String input, rows, output;
	private int quota=-1;
	private long seed=-1;
}
