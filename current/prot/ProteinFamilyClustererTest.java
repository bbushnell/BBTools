package prot;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Synthetic test for {@link ProteinFamilyClusterer}: two well-separated protein
 * families (hydrophobic vs charged), each with mutated members. Verifies the
 * 4-step pipeline recovers exactly two pure families, is order-independent (the
 * property the reassignment step buys), and produces a consensus that each member
 * matches better than the other family's consensus.
 *
 * <p>Run: {@code java -cp current --add-modules jdk.incubator.vector prot.ProteinFamilyClustererTest}</p>
 *
 * @author Eru
 */
public final class ProteinFamilyClustererTest {

	private static int passed=0, failed=0;

	//Family A: hydrophobic ILV backbone with scattered substitutions.
	private static final String[] FAM_A={
		"ILVILVILVILVILVILVILVILV",
		"ILVMLVILVMLVILVILVMLVILV",
		"MLVILVILVILAMLVILVILVILV",
		"ILVILAILVILVILVMLVILVILA",
	};
	//Family B: charged DEKR backbone with scattered substitutions.
	private static final String[] FAM_B={
		"DEKRDEKRDEKRDEKRDEKRDEKR",
		"DEKRDEQRDEKRDEKRDENRDEKR",
		"DENRDEKRDEKRDEKRDEKRDEQR",
		"DEKRDEKRDEQRDEKRDEKRDEKR",
	};

	public static void main(String[] args){
		final List<ProteinSequence> order1=build(new int[]{0,1,2,3, 4,5,6,7});//A then B
		final List<ProteinSequence> order2=build(new int[]{4,0,6,2, 5,3,7,1});//interleaved

		final List<ProteinFamily> f1=newClusterer().cluster(order1);
		final List<ProteinFamily> f2=newClusterer().cluster(order2);

		check("exactly two families recovered", f1.size()==2);
		check("each family is pure (single true family)", allPure(f1));
		check("all 8 sequences assigned", totalMembers(f1)==8);
		check("order-independent partition", samePartition(f1, f2));
		check("consensus is family-specific (members match own consensus best)",
				consensusIsFamilySpecific(f1));

		System.err.println("\n"+passed+" passed, "+failed+" failed.");
		if(failed>0){System.exit(1);}
	}

	private static ProteinFamilyClusterer newClusterer(){
		final ProteinFamilyClusterer c=new ProteinFamilyClusterer();
		c.clusterIdentity=30.0;
		c.reassignIdentity=30.0;
		c.minCoverage=0.6;
		c.verbose=false;
		return c;
	}

	/** Builds the 8-sequence input in the given permutation of indices 0-7. */
	private static List<ProteinSequence> build(int[] order){
		final String[] all=new String[8];
		System.arraycopy(FAM_A, 0, all, 0, 4);
		System.arraycopy(FAM_B, 0, all, 4, 4);
		final ArrayList<ProteinSequence> list=new ArrayList<ProteinSequence>();
		for(int idx : order){
			final String tag=(idx<4 ? "A"+idx : "B"+(idx-4));
			list.add(new ProteinSequence("p"+tag, all[idx]));
		}
		return list;
	}

	/** Every family's members must share the same true-family prefix letter (A/B). */
	private static boolean allPure(List<ProteinFamily> fams){
		for(ProteinFamily f : fams){
			char fam=0;
			for(ProteinSequence m : f.members){
				final char c=m.id.charAt(1);//"pA0" -> 'A'
				if(fam==0){fam=c;}
				else if(fam!=c){return false;}
			}
		}
		return true;
	}

	private static int totalMembers(List<ProteinFamily> fams){
		int n=0;
		for(ProteinFamily f : fams){n+=f.size();}
		return n;
	}

	/** Two clusterings agree iff their sets of member-id sets are identical. */
	private static boolean samePartition(List<ProteinFamily> a, List<ProteinFamily> b){
		return partitionSet(a).equals(partitionSet(b));
	}

	private static HashSet<String> partitionSet(List<ProteinFamily> fams){
		final HashSet<String> parts=new HashSet<String>();
		for(ProteinFamily f : fams){
			final ArrayList<String> ids=new ArrayList<String>();
			for(ProteinSequence m : f.members){ids.add(m.id);}
			java.util.Collections.sort(ids);
			parts.add(ids.toString());
		}
		return parts;
	}

	/** Each member must align to its own family's consensus at least as well as the other's. */
	private static boolean consensusIsFamilySpecific(List<ProteinFamily> fams){
		if(fams.size()!=2){return false;}
		final ProteinFamily x=fams.get(0), y=fams.get(1);
		return membersPreferOwn(x, y) && membersPreferOwn(y, x);
	}

	private static boolean membersPreferOwn(ProteinFamily own, ProteinFamily other){
		for(ProteinSequence m : own.members){
			final double self=pident(m, own.consensus);
			final double cross=pident(m, other.consensus);
			if(self<=cross){return false;}
		}
		return true;
	}

	private static double pident(ProteinSequence q, ProteinSequence t){
		final AAAlignment aln=AAAligner.align(q.enc, t.enc);
		return aln==null ? 0 : aln.pident();
	}

	private static void check(String name, boolean ok){
		System.err.println((ok ? "PASS" : "FAIL")+": "+name);
		if(ok){passed++;}else{failed++;}
	}
}
