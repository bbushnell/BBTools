package prok;

import java.util.Arrays;
import java.util.List;

import dna.AminoAcid;
import stream.Read;

/** One-genome, single-use CDS-inner adaptation. No calling, HBM loading, or global configuration.
 * The caller supplies every contig once and labels the returned final DP path. By default labels
 * are actual supported-stop evidence (not bonuses); the explicit hbmtrainingset=selected experiment
 * labels all selected CDS positive, without asserting HBM support. Completeness/truncation gating
 * remains inside this trainer. Biological boundary flags/code belong to the caller.
 * Whole selected complete CDSs are training hypotheses; stop-tail evidence does not validate starts.
 * The base must remain immutable until finish; no contig, ORF, or label array is retained.
 * @author Sayu
 */
public final class ProteinGuidedModel {

	public ProteinGuidedModel(final GeneModel base_){
		if(base_==null || base_.statsCDS.inner==null || base_.basesProcessed<0){
			throw new IllegalArgumentException("Inner adaptation requires a base model with nonnegative base count");
		}
		final FrameStats inner=base_.statsCDS.inner;
		if(inner.frames!=3 || inner.k<1 || inner.k>15 || inner.leftOffset!=0){
			throw new IllegalArgumentException("CDS adaptation requires three phases and ordinary inner-kmer geometry");
		}
		base=base_;
		delta=new FrameStats(inner.name, inner.k, inner.frames, inner.leftOffset);
	}

	/** Labels are parallel to finalSelected: normally evidence.stop.supported(), including zero-bonus
	 * support; hbmtrainingset=selected instead supplies true for every selected CDS. A true label is
	 * eligibility for positive training, not proof of support: incomplete/truncated CDS remain excluded.
	 * Coordinates must be restored genomic coordinates from GeneCaller.callGenes, not candidate lists.
	 * RNA entries have false labels and do not contribute CDS masks. Failed trainers cannot be reused. */
	public void addContig(final Read contig, final List<Orf> finalSelected, final boolean[] stopSupported){
		requireOpen();
		try{
			if(contig==null || contig.bases==null || contig.id==null || finalSelected==null
				|| stopSupported==null || finalSelected.size()!=stopSupported.length){
				throw new IllegalArgumentException("Contig and final-path ORFs require an equally sized supported-stop label vector");
			}
			final byte[] bases=contig.bases;
			for(byte b : bases){
				if(b<0){throw new IllegalArgumentException("Expected ASCII nucleotide bytes for inner adaptation");}
			}
			for(int j=0; j<finalSelected.size(); j++){
				final Orf orf=finalSelected.get(j);
				if(orf==null){throw new IllegalArgumentException("Null final-path feature");}
				if(orf.type!=ProkObject.CDS){
					if(stopSupported[j]){throw new IllegalArgumentException("Only CDS features can carry supported-stop labels");}
					continue;
				}
				if(orf.strand<0 || orf.strand>1 || orf.flipped()!=orf.strand || orf.scaflen!=bases.length
					|| !contig.id.equals(orf.scafName) || orf.start<0 || orf.stop<orf.start || orf.stop>=bases.length
					|| orf.frame<0 || orf.frame>2){
					throw new IllegalArgumentException("Final CDS must belong to this contig in restored genomic coordinates");
				}
				final int first=orf.strand==0 ? orf.start : bases.length-1-orf.stop;
				if(first%3!=orf.frame){throw new IllegalArgumentException("CDS phase disagrees with its oriented start");}
			}
			inputBases=Math.addExact(inputBases, bases.length);// Once per contig, not per strand or selected CDS.
			final byte[] positives=new byte[bases.length], excluded=new byte[bases.length];
			for(int strand=0; strand<2; strand++){
				Arrays.fill(positives, (byte)0); Arrays.fill(excluded, (byte)0);
				for(int j=0; j<finalSelected.size(); j++){
					final Orf orf=finalSelected.get(j);
					if(orf.type!=ProkObject.CDS || orf.strand!=strand){continue;}
					final int first=strand==0 ? orf.start : bases.length-1-orf.stop;
					final int last=strand==0 ? orf.stop : bases.length-1-orf.start;
					final boolean positive=stopSupported[j] && !orf.startTruncated && !orf.stopTruncated
						&& orf.length()%3==0;
					// A synthetic terminal boundary has no biological stop codon to trim off.
					final int end=orf.stopTruncated ? last : last-3;
					final byte[] mask=positive ? positives : excluded;
					// Same endpoint-relative phase convention as GeneModel.markFrames.
					for(long pos=(long)first+delta.k-1; pos<=end; pos++){
						final int i=(int)pos;
						mask[i]|=(byte)(1<<((i-first)%3));
					}
				}
				count(strand==0 ? bases : AminoAcid.reverseComplementBases(bases), positives, excluded);
			}
		}catch(RuntimeException | Error failure){consumed=true; throw failure;}
	}

	/** Count each unambiguous window once per phase. A positive bit takes precedence over exclusion. */
	private void count(final byte[] bases, final byte[] positives, final byte[] excluded){
		assert(bases.length==positives.length && bases.length==excluded.length)
			: "Endpoint masks must cover exactly the oriented contig being counted";
		int kmer=0, length=0;
		for(int i=0; i<bases.length; i++){
			final int x=AminoAcid.baseToNumber[bases[i]];
			if(x<0){length=0; kmer=0; continue;}
			kmer=((kmer<<2)|x)&delta.mask;
			if(++length<delta.k){continue;}
			for(int phase=0; phase<3; phase++){
				final int bit=1<<phase, valid=(positives[i]&bit)==0 ? 0 : 1;
				if(valid==0 && (excluded[i]&bit)!=0){continue;}
				delta.counts[valid][phase][kmer]=Math.addExact(delta.counts[valid][phase][kmer], 1);
				delta.validSums[valid]=Math.addExact(delta.validSums[valid], 1);
			}
		}
	}

	/** Destructively merge raw counts from a joined worker using the identical immutable base.
	 * Neither trainer may be accessed concurrently. A valid donor is consumed before any addition;
	 * an unrelated donor rejected by preconditions is untouched, but every failure consumes this receiver. No partial result may
	 * be reused. Merge even empty/background-only workers, then finish the receiver exactly once:
	 * weighting or rounding each worker separately would change the genome-wide training result. */
	public void merge(final ProteinGuidedModel donor){
		requireOpen();
		try{
			if(donor==null || donor==this){throw new IllegalArgumentException("Raw merge requires a distinct nonnull trainer");}
			donor.requireOpen();
			if(base!=donor.base || !delta.compatibleWith(donor.delta)){
				throw new IllegalArgumentException("Raw trainers must share the identical base model and inner geometry");
			}
			donor.consumed=true;
			if(inputBases<0 || donor.inputBases<0){throw new IllegalArgumentException("Negative raw input base count");}
			inputBases=Math.addExact(inputBases, donor.inputBases);
			for(int valid=0; valid<2; valid++){
				for(int phase=0; phase<delta.frames; phase++){
					final long[] target=delta.counts[valid][phase], source=donor.delta.counts[valid][phase];
					for(int kmer=0; kmer<delta.kMax; kmer++){
						if(target[kmer]<0 || source[kmer]<0){throw new IllegalArgumentException("Negative raw inner count");}
						target[kmer]=Math.addExact(target[kmer], source[kmer]);
					}
				}
				if(delta.validSums[valid]<0 || donor.delta.validSums[valid]<0){
					throw new IllegalArgumentException("Negative raw inner class total");
				}
				delta.validSums[valid]=Math.addExact(delta.validSums[valid], donor.delta.validSums[valid]);
			}
		}catch(RuntimeException | Error failure){consumed=true; throw failure;}
	}

	/** Base + weighted inner counts only. No usable positive windows means exact base identity.
	 * Unlike legacy FrameStats.multiplyBy, totals of the NEW delta are summed after cell rounding.
	 * Existing base totals are preserved, even when prior weighting made them differ from cell sums. */
	public GeneModel finish(){
		requireOpen(); consumed=true;
		if(delta.validSums[1]==0){return base;}
		final double weight=0.50*Math.max(100L, base.basesProcessed)/Math.max(100L, inputBases);
		if(!Double.isFinite(weight) || weight<0){throw new IllegalArgumentException("Nonfinite/negative inner adaptation weight");}
		final GeneModel out=base.copy();
		final FrameStats inner=out.statsCDS.inner;
		for(int valid=0; valid<2; valid++){
			long total=0;
			for(int phase=0; phase<3; phase++){
				for(int kmer=0; kmer<delta.kMax; kmer++){
					final long old=inner.counts[valid][phase][kmer], observed=delta.counts[valid][phase][kmer];
					final double scaled=observed*weight;
					if(old<0 || observed<0 || !Double.isFinite(scaled) || scaled>=0x1p63){
						throw new IllegalArgumentException("Inner counts must remain nonnegative representable longs");
					}
					final long addition=Math.round(scaled);
					inner.counts[valid][phase][kmer]=Math.addExact(old, addition);
					total=Math.addExact(total, addition);
				}
			}
			if(inner.validSums[valid]<0){throw new IllegalArgumentException("Negative base inner class total");}
			inner.validSums[valid]=Math.addExact(inner.validSums[valid], total);
		}
		// calculate() adds long counts before converting to double; guard both aggregate and cell sums.
		Math.addExact(inner.validSums[0], inner.validSums[1]);
		for(int phase=0; phase<3; phase++){
			for(int kmer=0; kmer<inner.kMax; kmer++){
				Math.addExact(inner.countsFalse[phase][kmer], inner.countsTrue[phase][kmer]);
			}
		}
		inner.calculate();
		return out;
	}

	private void requireOpen(){
		if(consumed){throw new IllegalStateException("Inner trainer was consumed by finish, donation, or failure");}
	}
	private final GeneModel base;
	private final FrameStats delta;
	private long inputBases;
	private boolean consumed;
}
