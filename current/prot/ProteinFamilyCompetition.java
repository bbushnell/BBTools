package prot;

import java.util.Arrays;
import java.util.HashSet;

/** Competitive assignment across independently stored, compatible family packs.
 * The caller supplies canonical pack order, preserving each pack's local rank.
 * One global F4 shortlist precedes positional alignment; the result identifies
 * its pack and local rank so labels retain their original package bindings.
 * Libraries are immutable; each caller thread owns a Worker.
 * @author Keqing
 */
public final class ProteinFamilyCompetition {

	/** Binds already verified packages. Rejects duplicate identities or different
	 * scoring backgrounds rather than comparing scores from incompatible spaces.
	 * This layer does not discover packs or decide which development packs to enable.
	 */
	public ProteinFamilyCompetition(final ProteinFamilyAssigner... packs_){
		require(packs_!=null && packs_.length>0, "Combined assignment requires at least one family pack");
		packs=packs_.clone();
		final HashSet<String> representatives=new HashSet<String>();
		long count=0;
		for(int p=0; p<packs.length; p++){
			require(packs[p]!=null, "Null family pack at index "+p);
			if(p>0){compatible(packs[0], packs[p], p);}
			count+=packs[p].familyCount();
			require(count<=Integer.MAX_VALUE, "Combined family count exceeds the supported index range");
			for(int r=0; r<packs[p].familyCount(); r++){
				require(representatives.add(packs[p].representative(r)),
					"Duplicate representative across packs: "+packs[p].representative(r));
			}
		}
		familyCount=(int)count;
		//One-time primitive sort avoids a boxed set and detects IDs shared by different representatives.
		final int[] ids=new int[familyCount]; int at=0;
		for(ProteinFamilyAssigner pack : packs){
			for(int r=0; r<pack.familyCount(); r++){ids[at++]=pack.familyId(r);}
		}
		Arrays.sort(ids);
		for(int i=1; i<ids.length; i++){require(ids[i]!=ids[i-1], "Duplicate permanent family ID across packs: "+ids[i]);}
		assert(familyCount>=SHORTLIST) : "Verified packages must supply enough families for one full top50";
	}

	public int familyCount(){return familyCount;}
	public int packCount(){return packs.length;}
	public Worker newWorker(){return new Worker();}

	/** Compares actual query transforms; profile beta/clipping/gap are already
	 * fixed and checked by each ProteinFamilyAssigner loader's core metadata.
	 */
	private static void compatible(final ProteinFamilyAssigner a, final ProteinFamilyAssigner b, final int p){
		final FamilyShortlistSidecar x=a.library().sidecar, y=b.library().sidecar;
		require(a.backgroundPin().equals(b.backgroundPin()), "Different positional background in pack "+p);
		require(x.dims==y.dims && x.compositionK==y.compositionK && x.coveringK==y.coveringK
			&& x.compositionAlphabet.groups().equals(y.compositionAlphabet.groups())
			&& x.coveringAlphabet.delegate.groups().equals(y.coveringAlphabet.delegate.groups())
			&& x.variant.equals(y.variant) && Double.compare(x.rawPseudocount, y.rawPseudocount)==0
			&& Arrays.equals(x.background, y.background),
			"Different composition or covering-kmer scoring recipe in pack "+p);
	}

	/** Reuses all shortlist/selection arrays; query results are detached snapshots. */
	public final class Worker {
		private Worker(){
			final int n=packs.length;
			shortlists=new ProteinSearcher.ShortlistScratch[n];
			selections=new HbmCompetitiveGate.Scratch[n];
			hits=new int[n][]; composition=new float[n][]; scores=new float[n][];
			medians=new double[n][]; spreads=new double[n][];
			localRanks=new int[n][SHORTLIST]; counts=new int[n]; winners=new int[n];
			for(int p=0; p<n; p++){
				final FamilyShortlistSidecar sidecar=packs[p].library().sidecar;
				final ProteinSearcher.ShortlistScratch s=new ProteinSearcher.ShortlistScratch(sidecar.dims, sidecar.nFamilies, SHORTLIST);
				shortlists[p]=s; hits[p]=s.sk; composition[p]=s.compScores; scores[p]=s.f4;
				medians[p]=sidecar.medLen; spreads[p]=sidecar.sdLnLen;
				selections[p]=new HbmCompetitiveGate.Scratch(sidecar.nFamilies, SHORTLIST);
			}
		}

		/** Scores only the global top50, then chooses the highest eligible positional
		 * score with the unchanged ASCII representative tie break. Counters cover all
		 * 50 candidates, even when the winner is in the first pack.
		 */
		public ProteinFamilyAssigner.Hit assign(final String queryId, final byte[] aminoAcids){
			final ProteinSequence query=new ProteinSequence(queryId, aminoAcids);
			for(int p=0; p<packs.length; p++){
				//The raw arrays are reused; local F4 is overwritten using global moments below.
				ProteinSearcher.scoreFamiliesF4(query, packs[p].library().sidecar, shortlists[p]);
			}
			FamilyShortlistScores.score(hits, composition, medians, spreads, query.length(), scores);
			FamilyShortlistScores.top(scores, topPacks, topRanks);
			Arrays.fill(counts, 0);
			for(int i=0; i<SHORTLIST; i++){final int p=topPacks[i]; localRanks[p][counts[p]++]=topRanks[i];}
			int best=-1, scored=0, rejected=0;
			for(int p=0; p<packs.length; p++){
				final HbmCompetitiveAssign.Library library=packs[p].library();
				final HbmCompetitiveGate.Scratch work=selections[p];
				final int rank=library.gate.select(query.enc, localRanks[p], counts[p], work);
				winners[p]=rank; scored+=work.scored; rejected+=work.lengthRejected;
				if(rank>=0 && (best<0 || HbmCompetitiveGate.better(work.score64, library.ids[rank],
						selections[best].score64, packs[best].library().ids[winners[best]]))){best=p;}
			}
			assert(scored+rejected==SHORTLIST) : "Combined gate must account for exactly the global top50";
			final int rank=best<0 ? -1 : winners[best];
			final HbmCompetitiveGate.Scratch work=selections[best<0 ? 0 : best];
			return new ProteinFamilyAssigner.Hit(best, rank, best<0 ? -1 : packs[best].familyId(rank),
				best<0 ? null : packs[best].representative(rank), query.length(), work, scored, rejected);
		}

		private final ProteinSearcher.ShortlistScratch[] shortlists;
		private final HbmCompetitiveGate.Scratch[] selections;
		private final int[][] hits, localRanks;
		private final float[][] composition, scores;
		private final double[][] medians, spreads;
		private final int[] counts, winners;
		private final int[] topPacks=new int[SHORTLIST], topRanks=new int[SHORTLIST];
	}

	private static void require(final boolean ok, final String reason){if(!ok){throw new IllegalArgumentException(reason);}}
	private static final int SHORTLIST=50;
	private final ProteinFamilyAssigner[] packs;
	private final int familyCount;
}
