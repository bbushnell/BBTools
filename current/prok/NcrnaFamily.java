package prok;

import consensus.BaseGraph;
import map.LongHashSet;

/**
 * One conserved-ncRNA family's loaded resources + per-family scavenger config,
 * as a shared (read-only after loading), family-agnostic bundle. Populated by
 * CallGenes.loadNcrnaResources(); consumed by GeneCaller, which builds one
 * NcrnaScavenger per bundle per instance (mirroring TrnaCaller's per-instance
 * construction -- NcrnaScavenger's TrnaKmerIndex is per-thread-mutable, so the
 * bundle itself is shared but each GeneCaller instance gets its own scavenger).
 *
 * <p>Two bundles may share the same kmerSet/kLong (e.g. srp_small and
 * srp_large both seed from one combined srp_17mers.fa, per Noire, 2026-08-23:
 * "mixing 120bp and 300bp sequences in one consensus doesn't hold together"
 * -- so the CONSENSUS libraries are separate, but the conserved-kmer SEED set
 * is shared across both size classes).
 * @author Neptune
 */
public class NcrnaFamily {

	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_){
		this(name_, library_, models_, modelNames_, kmerSet_, kLong_, minLen_, windowPad_,
				7, 60, true, 11f, 0.48f, 0.072f, 12);
	}

	public NcrnaFamily(String name_, byte[][] library_, BaseGraph[] models_, String[] modelNames_,
			LongHashSet kmerSet_, int kLong_, int minLen_, int windowPad_,
			int indexK_, int indexTopN_, boolean adaptive_,
			float adaptFloor_, float adaptTopFrac_, float adaptQFrac_, int fixedMinHits_){
		name=name_;
		library=library_;
		models=models_;
		modelNames=modelNames_;
		kmerSet=kmerSet_;
		kLong=kLong_;
		minLen=minLen_;
		windowPad=windowPad_;
		indexK=indexK_;
		indexTopN=indexTopN_;
		adaptive=adaptive_;
		adaptFloor=adaptFloor_;
		adaptTopFrac=adaptTopFrac_;
		adaptQFrac=adaptQFrac_;
		fixedMinHits=fixedMinHits_;
	}

	public final String name;
	public final byte[][] library;
	public final BaseGraph[] models;
	public final String[] modelNames;
	public final LongHashSet kmerSet;
	public final int kLong;
	public final int minLen;
	public final int windowPad;

	public final int indexK;
	public final int indexTopN;
	public final boolean adaptive;
	public final float adaptFloor;
	public final float adaptTopFrac;
	public final float adaptQFrac;
	public final int fixedMinHits;
}
