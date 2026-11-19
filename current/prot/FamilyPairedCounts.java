package prot;

import structures.ByteBuilder;

/** Worker-local paired outcomes for OFF, multiply and add1 on the same query.
 * Oracle membership refers to the strict pre-scalar shortlist union, not truth.
 * @author Sayu
 */
public final class FamilyPairedCounts{

	public FamilyPairedCounts(int families){
		require(families>0,"Paired counts require a nonempty live family roster");
		arms=new FamilyPanelRunner.PanelCounts[3];
		for(int a=0;a<3;a++){arms[a]=new FamilyPanelRunner.PanelCounts(families);}
	}

	/** Copies primitive outcomes immediately; the scorer may reuse all supplied arrays. */
	public void observe(int length,String[] terminals,int[] winners,int oracle,int membership,int unionSize){
		require(length>0 && terminals!=null && terminals.length==3 && winners!=null && winners.length==3,
			"Every paired observation must contain the same nonempty query and all three arms");
		require(unionSize>=0 && unionSize<=150 && unionSize<=arms[0].winners.length,
			"The oracle union contains at most three strict top50 lists of distinct live families");
		require((unionSize==0 && oracle==-1 && membership==0) ||
			(unionSize>0 && oracle>=0 && oracle<arms[0].winners.length && membership>0 && membership<8),
			"A nonempty retrieval union has one live oracle winner present in at least one arm");
		// Validate the whole observation before modifying counters.
		for(int a=0;a<3;a++){
			final int stage=FamilyPanelRunner.stage(terminals[a]),family=winners[a];
			require(family>=-1 && family<arms[a].winners.length && (family>=0)==(stage>=1 && stage<=3),
				"Terminal stage and live winner must agree for every paired arm");
		}
		for(int a=0;a<3;a++){
			final FamilyPanelRunner.PanelCounts c=arms[a];
			c.queries++;c.residues+=length;c.stages[FamilyPanelRunner.stage(terminals[a])]++;
			if(winners[a]>=0){c.accepted++;c.winners[winners[a]]++;}
		}
		for(int a=1;a<3;a++){
			final int state=(winners[0]>=0 ? 1 : 0)|(winners[a]>=0 ? 2 : 0);
			transitions[a-1][state]++;
			if(state==3){if(winners[0]==winners[a]){sameWinner[a-1]++;}else{changedWinner[a-1]++;}}
		}
		oracleMasks[membership]++;unionPairs+=unionSize;
	}

	public void add(FamilyPairedCounts other){
		require(other!=null && other!=this && arms[0].winners.length==other.arms[0].winners.length,
			"Worker reductions must preserve the same live family roster and avoid self-merging");
		for(int a=0;a<3;a++){arms[a].add(other.arms[a]);}
		for(int a=0;a<2;a++){
			for(int s=0;s<4;s++){transitions[a][s]=Math.addExact(transitions[a][s],other.transitions[a][s]);}
			sameWinner[a]=Math.addExact(sameWinner[a],other.sameWinner[a]);
			changedWinner[a]=Math.addExact(changedWinner[a],other.changedWinner[a]);
		}
		for(int i=0;i<8;i++){oracleMasks[i]=Math.addExact(oracleMasks[i],other.oracleMasks[i]);}
		unionPairs=Math.addExact(unionPairs,other.unionPairs);
	}

	public void validate(String panel){
		for(int a=0;a<3;a++){arms[a].validate(panel+"/"+ARM_NAMES[a]);}
		final long n=arms[0].queries;long masks=0;
		for(long count:oracleMasks){require(count>=0,"Oracle counts cannot be negative");masks=Math.addExact(masks,count);}
		require(masks==n && unionPairs>=n-oracleMasks[0] && unionPairs<=150L*(n-oracleMasks[0]),
			"Oracle masks conserve all queries and empty unions contribute no alignment pairs");
		for(int a=1;a<3;a++){
			final long[] t=transitions[a-1];
			require(arms[a].queries==n && arms[a].residues==arms[0].residues && t[0]+t[1]+t[2]+t[3]==n &&
				t[1]+t[3]==arms[0].accepted && t[2]+t[3]==arms[a].accepted &&
				sameWinner[a-1]+changedWinner[a-1]==t[3],
				"Paired transitions must reconcile both marginal acceptance counts and shared winners");
		}
	}

	/** One comparison row per ON arm; availability is explicitly conditional on the union. */
	public ByteBuilder transitions(String panel){
		validate(panel);final ByteBuilder out=new ByteBuilder();
		for(int a=1;a<3;a++){
			long both=0,gained=0,lost=0,neither=0;
			for(int mask=1;mask<8;mask++){
				final boolean off=(mask&1)!=0,on=(mask&(1<<a))!=0;
				if(off && on){both+=oracleMasks[mask];}else if(on){gained+=oracleMasks[mask];}
				else if(off){lost+=oracleMasks[mask];}else{neither+=oracleMasks[mask];}
			}
			out.append(panel).tab().append(ARM_NAMES[a]).tab().append(arms[0].queries);
			for(long count:transitions[a-1]){out.tab().append(count);}
			out.tab().append(sameWinner[a-1]).tab().append(changedWinner[a-1])
				.tab().append(oracleMasks[0]).tab().append(both).tab().append(gained).tab().append(lost).tab().append(neither).nl();
		}
		return out;
	}
	public long queries(){return arms[0].queries;}
	public long accepted(int arm){
		require(arm>=0 && arm<3,"Arm index must name OFF, multiply or add1");
		return arms[arm].accepted;
	}
	public long oracleCount(int mask){
		require(mask>=0 && mask<8,"Oracle mask uses only the three arm bits");
		return oracleMasks[mask];
	}
	public long unionPairs(){return unionPairs;}
	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	final FamilyPanelRunner.PanelCounts[] arms;
	private final long[][] transitions=new long[2][4];
	private final long[] sameWinner=new long[2],changedWinner=new long[2],oracleMasks=new long[8];
	private long unionPairs;
	static final String[] ARM_NAMES={"off","multiply","add1"};
	public static final String TRANSITION_HEADER="panel\tarm\tqueries\tboth_reject\toff_only_accept"+
		"\ton_only_accept\tboth_accept\tsame_winner\tchanged_winner\tempty_union"+
		"\toracle_both\toracle_new_on\toracle_lost_on\toracle_neither\n";
}
