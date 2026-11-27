package test.aligner;

import aligner.CovarianceModel;
import aligner.CovarianceModelBandedCyk;
import aligner.CovarianceModelLengthBands;
import aligner.CovarianceModelPlacementBands;
import aligner.CovarianceModelStateCorridors;

/** Independent three-state pair grammar tests profiling without timing claims.
 * Literal cell counts distinguish allocated initialization from scored cells.
 * @author Brian, Raiden
 */
public final class CovarianceModelSearchProfileTest{

	public static void main(final String[] args){
		check(args.length==0, "This finite search-profile fixture takes no files or timing parameters");
		final CovarianceModel model=CovarianceModelCykTest.pair();
		model.nodeType[0]=CovarianceModel.ROOT;
		model.nodeType[1]=CovarianceModel.MATP;
		model.nodeType[2]=CovarianceModel.END;
		model.consensusRight[0]=3;
		model.consensusLeft[2]=1;
		model.consensusRight[2]=2;
		final byte[] sequence={'A', 'T'};
		final CovarianceModelPlacementBands placement=CovarianceModelPlacementBands.window(model, 2,
			new int[]{0, 1, 2}, new int[]{0, 1, 2}, 2);
		final CovarianceModelStateCorridors corridor=new CovarianceModelStateCorridors(model,
			new short[12], new boolean[]{true, true, true});
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 1000, 1000);
		final CovarianceModelBandedCyk.SearchProfile profile=new CovarianceModelBandedCyk.SearchProfile(true);
		for(int repeat=0; repeat<3; repeat++){
			final CovarianceModelBandedCyk.EndpointResult original=scorer.searchOverlapping(sequence, 2, 1, 2, placement);
			final CovarianceModelBandedCyk.EndpointResult measured=scorer.searchOverlapping(sequence, 2, 1, 2, placement, null, -1, profile);
			check(original.score==measured.score && original.targetFrom==measured.targetFrom
				&& original.targetTo==measured.targetTo && original.contacts==measured.contacts
				&& original.unrestrictedScore==measured.unrestrictedScore && original.cells==measured.cells,
				"Instrumentation must preserve the ordinary score, endpoints, contacts, ties and allocated domain");
			// Six triangular cells in S and MP plus three END cells. MP updates
			// only d=2 at j=2; S and END update all six and three respectively.
			check(profile.representedCells==15L*(repeat+1) && profile.updatedCells==10L*(repeat+1),
				"The two-base literal grammar allocates15cells but updates10; these are different work measures");
		}
		final long oldCells=profile.representedCells, oldUpdates=profile.updatedCells;
		final CovarianceModelBandedCyk.EndpointResult learned=scorer.searchOverlapping(sequence, 2, 1, 2, placement, corridor);
		final CovarianceModelBandedCyk.EndpointResult measured=scorer.searchOverlapping(sequence, 2, 1, 2, placement, corridor, -1, profile);
		check(learned.score==measured.score && learned.targetFrom==measured.targetFrom
			&& learned.targetTo==measured.targetTo && learned.contacts==measured.contacts,
			"The same learned-layout path must be preserved under profiling and recycled-deck reuse");
		check(profile.searches==4 && profile.representedCells-oldCells==3 && profile.updatedCells-oldUpdates==3,
			"Zero-residual anchors retain exactly the root, emitting pair and END cell");
		check(profile.layoutCpu>=0 && profile.fillCpu>=0 && profile.selectionCpu>=0 && profile.geometryCpu>=0,
			"CPU clocks are nonnegative observations; a tiny fixture makes no performance claim");
		check(profile.lastResult==measured && profile.kernel.updatedCells[CovarianceModel.S]==19
			&& profile.kernel.updatedCells[CovarianceModel.MP]==4 && profile.kernel.updatedCells[CovarianceModel.E]==10,
			"Three full pair domains and one anchored domain retain actual results and independently counted per-state updates");
		checkKernel(profile, 12, 48);
		fullGrammar(CovarianceModelCykTest.branch(), 24, 30, 10);
		fullGrammar(CovarianceModelCykTest.linear(), 34, 57, 0);
		System.err.println("CM_SEARCH_PROFILE_FIXTURE_PASS");
	}
	/** Two bases give six triangular cells per non-END state, three per END.
	 * A B row proposes1+1+2+1+2+3split trials; some have impossible child scores. */
	private static void fullGrammar(final CovarianceModel model, final long updates,
			final long initialized, final long splits){
		final CovarianceModelBandedCyk scorer=new CovarianceModelBandedCyk(model,
			new CovarianceModelLengthBands(model, 2), 10000, 10000);
		final byte[] sequence={'A', 'T'};
		final CovarianceModelBandedCyk.EndpointResult plain=scorer.searchOverlapping(sequence, 2, 1, 2, null);
		final CovarianceModelBandedCyk.SearchProfile profile=new CovarianceModelBandedCyk.SearchProfile(true);
		final CovarianceModelBandedCyk.EndpointResult measured=scorer.searchOverlapping(sequence, 2, 1, 2, null, null, -1, profile);
		check(plain.score==measured.score && plain.targetFrom==measured.targetFrom && plain.targetTo==measured.targetTo
			&& plain.contacts==measured.contacts && plain.unrestrictedScore==measured.unrestrictedScore
			&& plain.unrestrictedFrom==measured.unrestrictedFrom && plain.unrestrictedTo==measured.unrestrictedTo,
			"Profiling every ordinary state or bifurcation must retain both selected and unrestricted answers");
		check(profile.updatedCells==updates && profile.kernel.splitTrials==splits,
			"Hand-counted recurrence cells and proposed bifurcation splits use separate domains");
		checkKernel(profile, model.states(), initialized);
	}
	private static void checkKernel(final CovarianceModelBandedCyk.SearchProfile profile,
			final long expectedStates, final long expectedClear){
		long states=0, cleared=0, updates=0, phaseCpu=profile.kernel.setupCpu;
		for(int type=0; type<=CovarianceModel.B; type++){
			states+=profile.kernel.states[type];
			cleared+=profile.kernel.clearedSlots[type];
			updates+=profile.kernel.updatedCells[type];
			phaseCpu+=profile.kernel.acquireCpu[type]+profile.kernel.clearCpu[type]
				+profile.kernel.recurrenceCpu[type]+profile.kernel.releaseCpu[type];
			check(profile.kernel.states[type]==profile.kernel.allocations[type]+profile.kernel.reuses[type],
				"Every filled state acquires exactly one score/contact deck, including zero-cell states");
		}
		check(states==expectedStates && cleared==expectedClear && cleared==profile.representedCells
			&& updates==profile.updatedCells && phaseCpu<=profile.fillCpu,
			"State counts, initialized slots and nonoverlapping phases must conserve the enclosing real fill");
	}
	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
}
