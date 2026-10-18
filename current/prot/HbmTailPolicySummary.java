package prot;

import java.util.HashMap;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import structures.ByteBuilder;

/** Reports a fixed score-threshold grid from existing tail diagnostics.
 * Keeps all non-score counterfactual conditions unchanged; does not realign,
 * select a cutoff, or interpret synthetic decoys as biological truth.
 * @author Keqing
 */
public final class HbmTailPolicySummary{

	public static void main(String[] args){
		try{run(args);}catch(Throwable t){t.printStackTrace(System.err); System.exit(1);}
	}

	private static void run(String[] args){
		final HashMap<String,String> opt=HmmComparisonData.options(args);
		for(String key : opt.keySet()){require(key.equals("in") || key.equals("out"), "Unknown option: "+key);}
		final Counts[][][] counts=new Counts[SOURCES.length][BINS.length][VARIANTS.length];
		final HashMap<String,int[]> seen=new HashMap<String,int[]>();
		final LineParser1 lp=new LineParser1('\t');
		final ByteFile in=ByteFile.makeByteFile(HmmComparisonData.required(opt, "in"), false);
		long rows=0;
		try{
			final byte[] header=in.nextLine();
			require(header!=null, "Tail diagnostic is empty"); lp.set(header);
			require(lp.terms()==27 && lp.termEquals("query", 0) && lp.termEquals("query_source", 1)
				&& lp.termEquals("variant", 3) && lp.termEquals("length", 4) && lp.termEquals("status", 5)
				&& lp.termEquals("query_terminal_identical", 15) && lp.termEquals("shortlist_rank", 22)
				&& lp.termEquals("terminal20_score64", 23) && lp.termEquals("counterfactual_pass", 26),
				"Expected the 27-column HbmTailRobustness counterfactual schema");
			for(byte[] line=in.nextLine(); line!=null; line=in.nextLine()){
				lp.set(line); require(lp.terms()==27, "Malformed diagnostic row: "+(rows+2));
				final int source=index(lp, 1, SOURCES), variant=index(lp, 3, VARIANTS), length=lp.parseInt(4);
				require(length>0 && length<=500, "Tail length lies outside the frozen 500aa experiment");
				final String id=lp.parseString(0), key=source+":"+id;
				require(!id.isEmpty(), "Diagnostic query ID must be nonempty");
				int[] mask=seen.get(key);
				if(mask==null){mask=new int[1]; seen.put(key, mask);}
				final int bit=1<<variant;
				require((mask[0]&bit)==0, "Duplicate query/variant row: "+key+"/"+VARIANTS[variant]); mask[0]|=bit;
				final boolean current=status(lp)==HbmTailGate.SUPPORTED, candidate=lp.termEquals("true", 26);
				require(candidate || lp.termEquals("false", 26), "Counterfactual flag must be true or false");
				final boolean missing=lp.termEquals("NA", 23);
				final int score=missing ? Integer.MIN_VALUE : lp.parseInt(23), identities=lp.parseInt(15), rank=lp.parseInt(22);
				require(identities>=0 && identities<=20 && rank>0, "Invalid terminal identities or shortlist rank");
				require(!candidate || (!missing && score>0), "Base candidate already requires positive terminal score");
				final int bin=length<=100 ? 1 : length<=200 ? 2 : length<=300 ? 3 : 4;
				add(counts, source, 0, variant, current, candidate, score, identities, rank);
				add(counts, source, bin, variant, current, candidate, score, identities, rank); rows++;
			}
		}finally{require(!in.close(), "Tail diagnostic input failed to close");}
		require(rows>0, "No diagnostic rows found");
		for(Map.Entry<String,int[]> e : seen.entrySet()){
			require(e.getValue()[0]==(1<<VARIANTS.length)-1, "Query does not contain all nine variants: "+e.getKey());
		}
		final ByteStreamWriter out=HmmComparisonData.writer(HmmComparisonData.required(opt, "out"));
		out.println("source\tlength_bin\tvariant\tminimum_score64\tqueries\tcurrent_pass\tcandidate_pass\tcandidate_terminal_identities_le2\tcandidate_outside_top50");
		final ByteBuilder row=new ByteBuilder();
		for(int s=0; s<SOURCES.length; s++){
			for(int b=0; b<BINS.length; b++){
				for(int v=0; v<VARIANTS.length; v++){
					final Counts c=counts[s][b][v]; if(c==null){continue;}
					for(int t=0; t<CUTOFFS.length; t++){
						out.println(row.clear().append(SOURCES[s]).tab().append(BINS[b]).tab().append(VARIANTS[v])
							.tab().append(CUTOFFS[t]).tab().append(c.rows).tab().append(c.current).tab().append(c.accepted[t])
							.tab().append(c.lowIdentity[t]).tab().append(c.outsideTop[t]));
					}
				}
			}
		}
		HmmComparisonData.close(out);
		System.err.println("TAIL_POLICY_SUMMARY_PASS rows="+rows+" queries="+seen.size()+" thresholds="+CUTOFFS.length+" diagnostic_only=true");
	}

	private static void add(Counts[][][] counts, int source, int bin, int variant,
			boolean current, boolean candidate, int score, int identities, int rank){
		assert(source>=0 && bin>=0 && variant>=0) : "Summary counters require validated source/bin/variant indexes";
		Counts c=counts[source][bin][variant];
		if(c==null){c=new Counts(); counts[source][bin][variant]=c;}
		c.rows++; if(current){c.current++;}
		for(int t=0; t<CUTOFFS.length; t++){
			if(candidate && score>=CUTOFFS[t]){
				c.accepted[t]++; if(identities<=2){c.lowIdentity[t]++;} if(rank>50){c.outsideTop[t]++;}
			}
		}
	}
	private static int index(LineParser1 lp, int field, String[] values){
		for(int i=0; i<values.length; i++){if(lp.termEquals(values[i], field)){return i;}}
		throw new IllegalArgumentException("Unknown diagnostic category: "+lp.parseString(field));
	}
	private static int status(LineParser1 lp){
		for(int i=0; i<HbmTailGate.STATUS_COUNT; i++){if(lp.termEquals(HbmTailGate.statusName(i), 5)){return i;}}
		throw new IllegalArgumentException("Unknown tail acceptance status: "+lp.parseString(5));
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Counts{
		long rows, current;
		final long[] accepted=new long[CUTOFFS.length], lowIdentity=new long[CUTOFFS.length], outsideTop=new long[CUTOFFS.length];
	}
	// Zero preserves the base counterfactual's strict score>0 rule; later cutoffs are inclusive.
	private static final int[] CUTOFFS={0, 128, 256, 512, 1024, 1536, 2048, 3072, 4096};
	private static final String[] SOURCES={"consensus_self", "input"}, BINS={"all", "1-100", "101-200", "201-300", "301-500"};
	private static final String[] VARIANTS={"original", "replace10", "replace20", "replace40", "shuffle20", "trim10", "trim20", "append10", "append20"};
	private HbmTailPolicySummary(){}
}
