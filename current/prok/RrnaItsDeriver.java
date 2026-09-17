package prok;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import shared.Shared;
import shared.Tools;
import stream.Read;
import structures.ByteBuilder;

/** Pure helper deriving ITS1/ITS2/combined-ITS intervals from called 18S/5.8S/LSU anchor
 * records. No file I/O, no global state, no production wiring by itself -- CallGenes
 * calls this once per read and renders GFF/FASTA from its output.
 * MVP policy resolved by Ganyu 2026-09-17 (see
 * /mnt/c/codex-lbl/Qiqi/workspace/reports/callgenes_r58_lsu_its_source_contract_design_20260917.md
 * for the full unresolved-policy design this MVP narrows):
 *   - ITS1 = 18S-to-5.8S gap, ITS2 = 5.8S-to-LSU gap, combined ITS = 18S-to-LSU gap
 *     (the combined gap naturally includes any intervening 5.8S bases -- no separate
 *     composition step is needed since it is the SAME gap-between-anchors formula
 *     applied to the outer pair).
 *   - Combined ITS is attempted directly between 18S and LSU regardless of whether a
 *     5.8S anchor is present.
 *   - Multi-anchor cardinality: simple deterministic nearest-compatible-neighbor greedy
 *     matching, no anchor reuse WITHIN one product type's own pairing pass (ITS1, ITS2,
 *     and combined are three independent passes; the same 5.8S anchor may participate
 *     in both its ITS1 pairing and its ITS2 pairing, since those are different passes).
 *   - Max span: DEFAULT_MAX_SPAN, a generous MVP placeholder, not a measured value.
 * @author G11
 */
public class RrnaItsDeriver{

	/*--------------------------------------------------------------*/
	/*----------------            Anchor             ----------------*/
	/*--------------------------------------------------------------*/

	/** A called anchor feature: genomic 0-based inclusive coordinates, normalized so
	 * start<=stop regardless of input order. */
	public static class Anchor{

		public Anchor(String seqid_, int start_, int stop_, byte strand_, String family_, String sourceId_, int sourceOrdinal_, float score_){
			seqid=seqid_;
			start=Tools.min(start_, stop_);
			stop=Tools.max(start_, stop_);
			strand=strand_;
			family=family_;
			sourceId=sourceId_;
			sourceOrdinal=sourceOrdinal_;
			score=score_;
		}

		public final String seqid;
		public final int start, stop;
		/** 0=plus, 1=minus; matches the existing prok strand convention. */
		public final byte strand;
		/** Explicit semantic identity: S18, R58, or LSU. Never inferred from a model name. */
		public final String family;
		public final String sourceId;
		public final int sourceOrdinal;
		public final float score;

		/** Overlap-taint scratch state, mirroring CallITS.Flank exactly (ported algorithm,
		 * 2026-09-17). Callers must build a FRESH Anchor per read/group -- these fields are
		 * mutated in place by resolveGroup() and are not meant to be reused across calls. */
		final boolean[] taintedWith=new boolean[3];//indexed by famIndex(): S18=0,R58=1,LSU=2
		boolean anyTaint(){return taintedWith[0]||taintedWith[1]||taintedWith[2];}
		boolean selfAmbiguous=false;
		String ambigWith=null;
	}

	private static int famIndex(String fam){
		return fam.equals(S18) ? 0 : fam.equals(R58) ? 1 : fam.equals(LSU) ? 2 : -1;
	}

	/*--------------------------------------------------------------*/
	/*----------------         Rejection Reason      ----------------*/
	/*--------------------------------------------------------------*/

	public static enum Reason{ OK, CROSS_CONTIG, MIXED_STRAND, OVERLAP, ZERO_LENGTH, MAX_SPAN }

	/*--------------------------------------------------------------*/
	/*----------------            Derived             ----------------*/
	/*--------------------------------------------------------------*/

	/** One successfully-derived interval. */
	public static class Derived{
		public Derived(String product_, String seqid_, int start_, int stop_, byte strand_, byte[] sequence_, Anchor from_, Anchor to_){
			product=product_; seqid=seqid_; start=start_; stop=stop_; strand=strand_; sequence=sequence_; from=from_; to=to_;
		}
		/** "ITS1", "ITS2", or "ITS" (combined). */
		public final String product;
		public final String seqid;
		/** genomic 0-based inclusive. */
		public final int start, stop;
		public final byte strand;
		/** Transcript orientation: reverse-complemented once, exactly, for minus strand. */
		public final byte[] sequence;
		/** Upstream and downstream bounding anchors, in transcript order. */
		public final Anchor from, to;

		/** Renders one GFF row for this derived interval. Deliberately NOT an Orf/
		 * Orf.appendGff() (Ganyu, 2026-09-17, source-inspection finding): Orf's
		 * constructor reads start/stop codons and assumes a >=3bp CDS-shaped feature,
		 * and adding derived ITS intervals to an ORF list would contaminate DP/stats.
		 * This mirrors Orf.appendGff()'s column layout and comma-style attrs by hand
		 * instead, at arm's length from the ORF machinery. Column 3 is the literal
		 * "internal_transcribed_spacer_region" (not a ProkObject type -- ITS is not a
		 * registered ProkObject type and this MVP does not add one). */
		public ByteBuilder appendGff(ByteBuilder bb){
			bb.append(truncatedSeqid());
			bb.tab();
			bb.append("BBTools").append('\t');
			bb.append("internal_transcribed_spacer_region").append('\t');
			bb.append(start+1).append('\t');
			bb.append(stop+1).append('\t');
			bb.append('.').append('\t');
			bb.append(Shared.strandCodes2[strand]).append('\t');
			bb.append('0').append('\t');
			bb.append("product:").append(product).append(',');
			bb.append("from:").append(from.sourceId).append(',');
			bb.append("to:").append(to.sourceId).append(',');
			bb.append("len:").append(sequence.length);
			return bb;
		}

		/** Builds a FASTA Read for outits=, transcript-oriented (sequence[] is already
		 * reverse-complemented for minus-strand derivations -- never flipped again here). */
		public Read toRead(){
			return new Read(sequence, null,
				truncatedSeqid()+"\t"+product+"\t"+Shared.strandCodes2[strand]+"\t"+(start+1)+"-"+(stop+1)
				+"\tfrom:"+from.sourceId+"\tto:"+to.sourceId, 0, 0);
		}

		/** GFF column 1 stops at the first whitespace (standard GFF seqid convention); this
		 * mirrors that for the FASTA header too, so a consumer can join GFF ITS rows back to
		 * outits FASTA records by contig name. Bug found live on a real Dori canary run,
		 * 2026-09-17 (G11): appendGff() already truncated at print time but toRead() used the
		 * caller-supplied readId verbatim, so any genome whose FASTA headers carry descriptive
		 * text after the accession (i.e. most real genomes) produced a FASTA header that could
		 * never be joined back to its own GFF row by contig name. */
		private String truncatedSeqid(){
			for(int i=0, max=seqid.length(); i<max; i++){
				char c=seqid.charAt(i);
				if(c==' ' || c=='\t'){return seqid.substring(0, i);}
			}
			return seqid;
		}
	}

	/** Either a successful Derived record (reason==OK) or a rejection reason with no
	 * Derived -- rejected pairs never produce a GFF row or outits record. */
	public static class Result{
		Result(Derived derived_){ derived=derived_; reason=Reason.OK; }
		Result(Reason reason_){ derived=null; reason=reason_; assert(reason_!=Reason.OK) : "use the Derived constructor for OK"; }
		public final Derived derived;
		public final Reason reason;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Constants            ----------------*/
	/*--------------------------------------------------------------*/

	/** MVP placeholder (Ganyu, 2026-09-17): a generous span cap so the mechanism exists
	 * and is easy to re-parameterize once a real value is measured/approved. Not a
	 * scientifically derived default. */
	public static final int DEFAULT_MAX_SPAN=10000;

	public static final String S18="s18", R58="r58", LSU="lsu";

	/*--------------------------------------------------------------*/
	/*----------------         Pairwise Derive       ----------------*/
	/*--------------------------------------------------------------*/

	/** Attempts to derive the gap interval between two specific anchors. Determines
	 * transcript order internally from the strand convention (plus: lower coordinate is
	 * upstream; minus: higher coordinate is upstream) -- callers pass anchors in either
	 * order. Never guesses past an invalid geometry; returns a rejection Result instead. */
	public static Result derive(String product, Anchor a, Anchor b, byte[] bases, int maxSpan){
		if(!a.seqid.equals(b.seqid)) return new Result(Reason.CROSS_CONTIG);
		if(a.strand!=b.strand) return new Result(Reason.MIXED_STRAND);

		final int ovLo=Tools.max(a.start, b.start), ovHi=Tools.min(a.stop, b.stop);
		if(ovHi>=ovLo) return new Result(Reason.OVERLAP);

		//gap formula from the design doc, symmetric in a/b: gapStart=min(stops)+1, gapStop=max(starts)-1
		final int gapStart=Tools.min(a.stop, b.stop)+1;
		final int gapStop=Tools.max(a.start, b.start)-1;
		if(gapStart>gapStop) return new Result(Reason.ZERO_LENGTH);

		final int gapLen=gapStop-gapStart+1;
		if(maxSpan>=0 && gapLen>maxSpan) return new Result(Reason.MAX_SPAN);

		final Anchor upstream=(a.strand==0) ? (a.start<b.start ? a : b) : (a.start>b.start ? a : b);
		final Anchor downstream=(upstream==a) ? b : a;

		byte[] seq=Arrays.copyOfRange(bases, gapStart, gapStop+1);
		Tools.toUpperCase(seq);
		if(a.strand==1){ seq=AminoAcid.reverseComplementBases(seq); }

		return new Result(new Derived(product, a.seqid, gapStart, gapStop, a.strand, seq, upstream, downstream));
	}

	/*--------------------------------------------------------------*/
	/*------  Multi-Anchor Pairing (ported from CallITS, proven)  ---*/
	/*--------------------------------------------------------------*/
	// Algorithm ported faithfully from prok/CallITS.java's markOverlaps()+partner()
	// (Brian/G11, 2026-09-09; acceptance authority = Qiqi's
	// fungal_its_synthetic_acceptance_fixtures_20260909.md), translated from GFF-parsed
	// Flank records to live Anchor objects for in-process CallGenes wiring (Ganyu,
	// 2026-09-17: "reuse code where practical; do not broaden scope by refactoring all of
	// CallITS today"). Semantics preserved exactly: mutual-nearest-neighbor pairing with
	// explicit equal-distance tie detection in BOTH directions (a tie or a non-mutual
	// nearest = no pairing, never first-wins); same-family overlap = unresolved copy
	// identity (both anchors excluded from ANY pairing); cross-family overlap = taints
	// only the regions touching that family pair; a full (combined-ITS) span is rejected
	// if it bridges over any overlap-tainted anchor.

	/** Marks overlap taint on every pair of intersecting anchors in one contig+strand
	 * group (all families together), matching CallITS.markOverlaps() exactly. */
	private static void markOverlaps(ArrayList<Anchor> group){
		for(int i=0; i<group.size(); i++){
			for(int j=i+1; j<group.size(); j++){
				Anchor a=group.get(i), b=group.get(j);
				if(a.start<=b.stop && b.start<=a.stop){//intersecting
					if(a.family.equals(b.family)){
						a.selfAmbiguous=b.selfAmbiguous=true;
						a.ambigWith=(a.ambigWith==null ? b.sourceId : a.ambigWith+"|"+b.sourceId);
						b.ambigWith=(b.ambigWith==null ? a.sourceId : b.ambigWith+"|"+a.sourceId);
					}else{
						a.taintedWith[famIndex(b.family)]=true;
						b.taintedWith[famIndex(a.family)]=true;
					}
				}
			}
		}
	}

	/** Finds the mutual-nearest partner of family `fam` for source anchor `a`, on the
	 * given genomic side (rightSide=true means the candidate must be genomic-right of a).
	 * Returns null (no pairing) on: no candidate, a tie at the minimum distance, a's own
	 * unresolved same-family ambiguity, the best candidate itself being unresolved, or a
	 * failed mutual-nearest check (including a reverse-direction tie). Ported from
	 * CallITS.partner() exactly -- never first-wins on a tie. */
	private static Anchor partner(ArrayList<Anchor> group, Anchor a, String fam, boolean rightSide){
		if(a.selfAmbiguous) return null;
		Anchor best=null; long bd=-1; int ties=0;
		for(Anchor c : group){
			if(!c.family.equals(fam)) continue;
			final long d=(rightSide ? c.start-(long)a.stop : a.start-(long)c.stop);
			if(d<=0) continue;
			if(bd<0 || d<bd){bd=d; best=c; ties=1;}
			else if(d==bd){ties++;}
		}
		if(best==null || ties>1 || best.selfAmbiguous) return null;

		Anchor rbest=null; long rbd=-1; int rties=0;
		for(Anchor c : group){
			if(!c.family.equals(a.family)) continue;
			final long d=(rightSide ? best.start-(long)c.stop : c.start-(long)best.stop);
			if(d<=0) continue;
			if(rbd<0 || d<rbd){rbd=d; rbest=c; rties=1;}
			else if(d==rbd){rties++;}
		}
		if(rbest!=a || rties>1) return null;
		return best;
	}

	/** Resolves ITS1/ITS2/combined-ITS for ALL anchors in one (seqid,strand) group
	 * (mixed families together, as markOverlaps needs to see every family to diagnose
	 * cross-family taint). Ambiguous or rejected pairs simply produce no Derived record
	 * -- never a fabricated or first-wins gap. */
	public static ArrayList<Derived> resolveGroup(ArrayList<Anchor> group, byte[] bases, boolean emitFull, boolean emitIndividual, int maxSpan){
		ArrayList<Derived> out=new ArrayList<Derived>();
		if(group==null || group.isEmpty()) return out;
		markOverlaps(group);
		final boolean plus=(group.get(0).strand==0);

		if(emitFull){
			for(Anchor a : group){
				if(!a.family.equals(S18)) continue;
				Anchor lsu=partner(group, a, LSU, plus);//transcript-downstream side holds LSU
				if(lsu==null) continue;
				if(a.taintedWith[famIndex(LSU)] || lsu.taintedWith[famIndex(S18)]) continue;//OVERLAPPING_FLANK
				Result r=derive("ITS", a, lsu, bases, maxSpan);
				if(r.reason!=Reason.OK) continue;
				//no-bridge rule: reject if the span contains any overlap-tainted anchor
				boolean bridgesTaint=false;
				for(Anchor g2 : group){
					if(g2.anyTaint() && g2.start<=r.derived.stop && r.derived.start<=g2.stop){bridgesTaint=true; break;}
				}
				if(bridgesTaint) continue;
				//both-spacers-empty rule: a 5.8S sitting immediately adjacent to BOTH outer
				//flanks means the span contains only 5.8S, not a real ITS1+5.8S+ITS2 product
				boolean bothSpacersEmpty=false;
				for(Anchor g2 : group){
					if(g2.family.equals(R58) && g2.start>=r.derived.start && g2.stop<=r.derived.stop){
						final Anchor gl=(plus ? a : lsu), gr=(plus ? lsu : a);
						if(g2.start-gl.stop<=1 && gr.start-g2.stop<=1){bothSpacersEmpty=true; break;}
					}
				}
				if(bothSpacersEmpty) continue;
				out.add(r.derived);
			}
		}

		if(emitIndividual){
			for(Anchor a : group){
				if(!a.family.equals(R58)) continue;
				//ITS1: 18S transcript-upstream of the 5.8S pivot
				if(!a.taintedWith[famIndex(S18)]){
					Anchor s18=partner(group, a, S18, !plus);
					if(s18!=null && !s18.taintedWith[famIndex(R58)]){
						Result r=derive("ITS1", a, s18, bases, maxSpan);
						if(r.reason==Reason.OK) out.add(r.derived);
					}
				}
				//ITS2: LSU transcript-downstream of the 5.8S pivot
				if(!a.taintedWith[famIndex(LSU)]){
					Anchor lsu=partner(group, a, LSU, plus);
					if(lsu!=null && !lsu.taintedWith[famIndex(R58)]){
						Result r=derive("ITS2", a, lsu, bases, maxSpan);
						if(r.reason==Reason.OK) out.add(r.derived);
					}
				}
			}
		}
		return out;
	}
}
