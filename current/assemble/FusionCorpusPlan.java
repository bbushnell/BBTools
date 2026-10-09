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

/** Reproducibly randomizes depth, error quality and multi-K recipes per genome. @author Fischl */
final class FusionCorpusPlan {

	/** Reads the verified seven-column reference manifest and appends independent random draws. */
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
		final ByteStreamWriter output=new ByteStreamWriter(parser.out1, false, false, false);
		output.start();
		final LineParser1 fields=new LineParser1('\t');
		final ByteBuilder row=new ByteBuilder();
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
			inputError=input.close();
			outputError=output.poisonAndWait();
		}
		if(inputError || outputError || genomes!=expected){
			throw new IllegalStateException("Incomplete corpus plan: rows="+genomes+", expected="+expected);
		}
		writeDistribution(summary, depths, qualities, lengths);
		System.err.println("FUSION_CORPUS_PLAN_PASS genomes="+genomes+" seed="+seed);
	}

	/** Writes every possible bin, including zeros, so missing range coverage is visible. */
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
