package assemble;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import parse.Parser;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Randomizes depth, error quality and multi-K recipes for a fixed ordered manifest.
 * One seeded random stream supplies all draws; changing row order or membership changes
 * later assignments. This planner preserves supplied metadata without verifying genome
 * contents, biological groups or representativeness, and does not simulate reads itself.
 * @author Fischl
 */
final class FusionCorpusPlan{

	/**
	 * Called by FusionJoinDiagnostic plan mode after PreParser/config processing. Skips blank
	 * and comment lines, requires seven columns and literal-unique IDs/source paths, and keeps
	 * each original row intact. Other metadata and the referenced genome files are not validated.
	 * Depth is drawn from 15–150, quality from 25–35. Sorted schedules contain one K in 24–48,
	 * one or two distinct Ks in 49–96, and one in 97–127; read seeds are nonnegative longs.
	 * The exact expected row count is checked after closing the plan, before the distribution
	 * file is written. A failure can leave a partial plan; the two outputs are not atomic.
	 * Initialization failures close acquired streams and preserve the original unchecked
	 * failure, with cleanup failures suppressed; partial output files are not removed.
	 * Processing cleanup always attempts writer shutdown. An unchecked input-close failure
	 * retains precedence over processing errors, with writer-cleanup failures suppressed.
	 * @param args Native input/output options, optional seed and expected row count
	 */
	static void run(final String[] args){
		final Parser parser=new Parser();
		long seed=350194;
		int expected=1000;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			String a=split[0].toLowerCase(java.util.Locale.ROOT);
			while(a.startsWith("-")){a=a.substring(1);}
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("seed")){seed=Long.parseLong(b);}
			else if(a.equals("expected")){expected=Integer.parseInt(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown corpus-plan option: "+arg);}
		}
		if(parser.in1==null || parser.out1==null || expected<1){
			throw new IllegalArgumentException("Plan mode requires in=genomes.tsv out=samples.tsv expected=1000 seed=N.");
		}
		final String summary=parser.out1+".distribution.tsv";
		if(!Tools.testInputFiles(false, true, parser.in1) ||
				!Tools.testOutputFiles(false, false, false, parser.out1, summary)){
			throw new IllegalArgumentException("Corpus input must exist and outputs must be new writable files.");
		}
		TadpoleGraph.checkPaths(parser.in1, new java.util.ArrayList<String>(), parser.out1, summary);
		final Random random=new Random(seed);
		final int[] depths=new int[151], qualities=new int[36], lengths=new int[128];
		final HashSet<String> ids=new HashSet<String>(), paths=new HashSet<String>();
		final ByteFile input=ByteFile.makeByteFile(parser.in1, true);
		ByteStreamWriter acquired=null;
		final LineParser1 fields;
		final ByteBuilder row;
		try{
			acquired=new ByteStreamWriter(parser.out1, false, false, false);
			acquired.start();
			fields=new LineParser1('\t');
			row=new ByteBuilder();
		}catch(RuntimeException|Error failure){
			closeAfterInitializationFailure(input, acquired, failure);
			throw failure;
		}
		final ByteStreamWriter output=acquired;
		int genomes=0;
		boolean inputError, outputError;
		try{
			output.print("#genome\tsource\ttaxid\tgroup\tcontig_length50\tbases\tscaffolds\tdepth\tquality\tk_schedule\tread_seed\tplan_seed\n");
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				if(line.length==0 || line[0]=='#'){continue;}
				fields.set(line);
				if(fields.terms()!=7 || !ids.add(fields.parseString(0)) || !paths.add(fields.parseString(1))){
					throw new IllegalArgumentException("Corpus references must have seven columns and unique IDs/paths.");
				}
				final int depth=15+random.nextInt(136), quality=25+random.nextInt(11);
				final int low=24+random.nextInt(25), middle=49+random.nextInt(48), high=97+random.nextInt(31);
				final int[] ks;
				if(random.nextBoolean()){
					int extra=49+random.nextInt(47);
					if(extra>=middle){extra++;}
					ks=new int[]{low, middle, extra, high};
				}else{ks=new int[]{low, middle, high};}
				Arrays.sort(ks);
				assert(ks[0]>=24 && ks[ks.length-1]<=127 && high<150) :
						"All corpus Ks must fit the150-base reads and preserve a long-read overlap margin.";
				row.clear().append(line).tab().append(depth).tab().append(quality).tab();
				for(int i=0; i<ks.length; i++){
					assert(i==0 || ks[i]>ks[i-1]) : "Random multi-K schedules must contain distinct lengths.";
					if(i>0){row.comma();}
					row.append(ks[i]);
					lengths[ks[i]]++;
				}
				row.tab().append(random.nextLong()&Long.MAX_VALUE).tab().append(seed).nl();
				output.print(row);
				depths[depth]++;
				qualities[quality]++;
				genomes++;
			}
		}finally{
			try{inputError=input.close();}
			catch(RuntimeException|Error failure){
				// ASM024: retain input-close precedence while still closing and joining the writer.
				try{
					if(output.poisonAndWait()){failure.addSuppressed(new IllegalStateException("Corpus plan output cleanup failed: "+output.fname));}
				}catch(RuntimeException|Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
				throw failure;
			}
			outputError=output.poisonAndWait();
		}
		if(inputError || outputError || genomes!=expected){
			throw new IllegalStateException("Incomplete corpus plan: rows="+genomes+", expected="+expected);
		}
		writeDistribution(summary, depths, qualities, lengths);
		System.err.println("FUSION_CORPUS_PLAN_PASS genomes="+genomes+" seed="+seed);
	}

	/** Closes resources acquired before the processing try; ASM024 output-open failure
	 * is reproduced by the public planner fixture. Writer cleanup is attempted even
	 * if input cleanup fails. State distinguishes unstarted ownership from a worker. */
	private static void closeAfterInitializationFailure(final ByteFile input, final ByteStreamWriter writer, final Throwable failure){
		assert(input!=null && failure!=null) : "run() calls initialization cleanup only after acquiring input and catching its primary failure.";
		try{
			if(input.close()){failure.addSuppressed(new IllegalStateException("Corpus manifest input cleanup failed."));}
		}catch(RuntimeException|Error cleanup){
			if(cleanup!=failure){failure.addSuppressed(cleanup);}
		}finally{
			if(writer!=null){
				try{
					final boolean failed=writer.getState()==Thread.State.NEW ? writer.closeBeforeStart() : writer.poisonAndWait();
					if(failed){failure.addSuppressed(new IllegalStateException("Corpus plan output cleanup failed: "+writer.fname));}
				}catch(RuntimeException|Error cleanup){
					if(cleanup!=failure){failure.addSuppressed(cleanup);}
				}
			}
		}
	}

	/**
	 * Writes all bins, including zeros. Depth/quality counts each sum to the number of genomes;
	 * K counts sum to the number of selected schedule entries, so each genome contributes three
	 * or four observations there. These are realized draw counts, not expected probabilities.
	 */
	private static void writeDistribution(final String path, final int[] depths,
			final int[] qualities, final int[] lengths){
		assert(depths.length==151 && qualities.length==36 && lengths.length==128) :
				"Distribution arrays must cover the full declared sampling ranges.";
		final ByteStreamWriter writer=new ByteStreamWriter(path, false, false, false);
		writer.start();
		boolean error;
		try{
			writer.print("parameter\tvalue\tcount\n");
			final ByteBuilder row=new ByteBuilder();
			for(int d=15; d<=150; d++){
				writer.print(row.clear().append("depth\t").append(d).tab().append(depths[d]).nl());
			}
			for(int q=25; q<=35; q++){
				writer.print(row.clear().append("quality\t").append(q).tab().append(qualities[q]).nl());
			}
			for(int k=24; k<=127; k++){
				writer.print(row.clear().append("k\t").append(k).tab().append(lengths[k]).nl());
			}
		}finally{error=writer.poisonAndWait();}
		if(error){throw new IllegalStateException("Could not save corpus distribution: "+path);}
	}
}
