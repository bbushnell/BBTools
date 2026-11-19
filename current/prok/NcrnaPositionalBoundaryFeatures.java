package prok;

import java.util.Arrays;

import consensus.BaseGraph;
import idaligner.AlignmentStats;

/** Worker-owned positional features shared by pooled-member training/inference.
 * Tables pool by reported member, while length and identity use the selected
 * candidate's actual consensus. Reuses the existing positional feature/Quantum
 * implementation unchanged; no truth coordinates participate in fill().
 * Scratch and result arrays are overwritten by the next call and never shared.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaPositionalBoundaryFeatures{

	/** Borrows finalized member tables, one for each end with identical radius. */
	public NcrnaPositionalBoundaryFeatures(final RrnaPositionalKmerTable.Table five,
			final RrnaPositionalKmerTable.Table three){
		this(five, three, null, false);
	}

	/** Optional three-probe layout, ordered best placement, left exclusive, right
	 * exclusive. Training owns LOO views; inference borrows only finalized scores.
	 */
	public NcrnaPositionalBoundaryFeatures(final RrnaPositionalKmerTable.Table[] five,
			final RrnaPositionalKmerTable.Table[] three, final boolean leaveOneOut){
		this(first(five), first(three), new RrnaPositionalKmerTable.Table[][]{five.clone(), three.clone()}, leaveOneOut);
	}

	private NcrnaPositionalBoundaryFeatures(final RrnaPositionalKmerTable.Table five,
			final RrnaPositionalKmerTable.Table three, final RrnaPositionalKmerTable.Table[][] probes_,
			final boolean leaveOneOut){
		this(five, three, probes_, null, leaveOneOut);
	}

	/** Full-error training/evaluation uses the same features with measured backgrounds. */
	public NcrnaPositionalBoundaryFeatures(final RrnaPositionalKmerTable.Table[] five,
			final RrnaPositionalKmerTable.Table[] three, final NcrnaWeightedEndpointTable[] fiveBackground,
			final NcrnaWeightedEndpointTable[] threeBackground, final boolean leaveOneOut){
		this(first(five), first(three), new RrnaPositionalKmerTable.Table[][]{five.clone(), three.clone()},
			new NcrnaWeightedEndpointTable[][]{fiveBackground.clone(), threeBackground.clone()}, leaveOneOut);
	}

	private NcrnaPositionalBoundaryFeatures(final RrnaPositionalKmerTable.Table five,
			final RrnaPositionalKmerTable.Table three, final RrnaPositionalKmerTable.Table[][] probes_,
			final NcrnaWeightedEndpointTable[][] backgrounds_, final boolean leaveOneOut){
		if(five==null || three==null || !five.finished || !three.finished
				|| !five.model.equals(three.model) || !five.end.equals("5prime") || !three.end.equals("3prime")
				|| five.radius!=three.radius){
			throw new IllegalArgumentException("Positional features require finalized same-member tables for both ends and one radius");
		}
		tables=new RrnaPositionalKmerTable.Table[]{five, three};
		probes=probes_;
		backgrounds=backgrounds_;
		radius=five.radius;
		final int width=probes==null ? five.inputs : Math.addExact(Math.multiplyExact(3, five.sites), 6);
		values=new float[][]{new float[width], new float[width]};
		fuzzStats=probes==null ? null : new AlignmentStats(true);
		loo=leaveOneOut ? new NcrnaPositionalLeaveOneOut[2][3] : null;
		if(probes!=null){
			for(int end=0; end<2; end++){
				if(backgrounds!=null && backgrounds[end].length!=3){
					throw new IllegalArgumentException("Each full-error end needs three matching weighted probes");
				}
				for(int probe=0; probe<3; probe++){
					final RrnaPositionalKmerTable.Table table=probes[end][probe];
					if(table==null || !table.finished || table.radius!=radius || table.k!=9
							|| !table.model.equals(tables[end].model) || !table.end.equals(tables[end].end)
							|| (probe==1 && table.offset!=-9) || (probe==2 && table.offset!=1)){
						throw new IllegalArgumentException("Three K9 probes must share member/end/radius and use explicit exclusive offsets -9 and +1");
					}
					if(backgrounds!=null){
						if(backgrounds[end][probe]==null || backgrounds[end][probe].geometry!=table){
							throw new IllegalArgumentException("Weighted probes must borrow these exact finalized correction tables");
						}
						backgrounds[end][probe].scoreKey(0);
					}
					if(loo!=null){
						loo[end][probe]=new NcrnaPositionalLeaveOneOut(table, backgrounds==null ? null : backgrounds[end][probe]);
					}
				}
			}
		}
	}

	private static RrnaPositionalKmerTable.Table first(final RrnaPositionalKmerTable.Table[] probes){
		if(probes==null || probes.length!=3){throw new IllegalArgumentException("Each end requires exactly three ordered positional probes");}
		return probes[0];
	}

	/** Bind a counted training locus once before its jittered training/monitor rows. */
	public void setLeaveOneOut(final byte[] bases, final int trueStart, final int trueStop){
		if(loo==null){throw new IllegalStateException("Only the optional training view excludes a known locus");}
		for(int end=0; end<2; end++){
			for(final NcrnaPositionalLeaveOneOut view:loo[end]){view.set(bases, end==0 ? trueStart : trueStop);}
		}
	}

	/** Fills candidate-span features in zero-based inclusive sense coordinates.
	 * @param fixedMask Bit1 fixes the5prime end; bit2 fixes the3prime end
	 * @return Usable-end bits in the same convention; a missing word keeps an end
	 * unavailable. Unusable arrays retain diagnostic poison and must never be read;
	 * only the returned bits establish availability, never a numeric feature value.
	 */
	public int fill(final byte[] bases, final int start, final int stop,
			final byte[] selectedConsensus, final float contigGC, final int fixedMask){
		if(probes!=null){throw new IllegalStateException("Three-probe features require the selected HBM and its coverage normalization");}
		if((fixedMask&~3)!=0){throw new IllegalArgumentException("Positional fixed mask contains only the two gene-end bits");}
		for(int end=0; end<2; end++){Arrays.fill(values[end], Float.NaN);}
		if(fixedMask==3){return 0;}
		validateGC(contigGC);
		final float identity=RrnaEndpointVector.candidateIdentity(bases, start, stop, selectedConsensus, stats);
		validateIdentity(identity);
		int usable=0;
		for(int end=0; end<2; end++){
			final int bit=1<<end;
			if((fixedMask&bit)!=0){continue;}
			final RrnaPositionalKmerTable.Table table=tables[end];
			if(!scoreWindow(table, bases, end==0 ? start : stop, values[end], 0)){continue;}
			values[end][table.sites]=(stop-start+1)/(float)selectedConsensus.length;
			values[end][table.sites+1]=contigGC;
			values[end][table.sites+2]=identity;
			validateValues(values[end]);
			usable|=bit;
		}
		assert((usable&fixedMask)==0) : "Physically fixed gene ends cannot enter positional inference or training rows";
		return usable;
	}

	/** Optional input order: three score windows, length ratio, contig GC, Quantum
	 * candidate identity, then sliding-defined HBM coverage at endpoint -1/0/+1.
	 * Both ends use this SAME unmoved span. Truth coordinates are confined to the
	 * separate training-only LOO binding and never select the fuzziness anchor.
	 */
	public int fill(final byte[] bases, final int start, final int stop, final byte[] selectedConsensus,
			final BaseGraph selectedModel, final float plateau, final float contigGC, final int fixedMask){
		if(probes==null){return fill(bases, start, stop, selectedConsensus, contigGC, fixedMask);}
		if((fixedMask&~3)!=0){throw new IllegalArgumentException("Positional fixed mask contains only the two gene-end bits");}
		validateGC(contigGC);
		for(int end=0; end<2; end++){Arrays.fill(values[end], Float.NaN);}
		if(fixedMask==3){return 0;}
		final float identity=RrnaEndpointVector.candidateIdentity(bases, start, stop, selectedConsensus, stats);
		validateIdentity(identity);
		final int sites=tables[0].sites, extras=3*sites;
		final byte[] candidate=Arrays.copyOfRange(bases, start, stop+1);
		TrnaBoundaryFeatures.tipFuzzinessFeatures(candidate, selectedModel, plateau, fuzzStats,
			values[0], extras+3, values[1], extras+3);
		for(int end=0; end<2; end++){
			for(int i=extras+3; i<extras+6; i++){
				if(!Float.isFinite(values[end][i])){throw new IllegalArgumentException("Computed HBM coverage features must be finite");}
			}
		}
		int usable=0;
		for(int end=0; end<2; end++){
			final int bit=1<<end;
			if((fixedMask&bit)!=0){
				Arrays.fill(values[end], Float.NaN);
				continue;
			}
			boolean complete=true;
			final int rawEnd=end==0 ? start : stop;
			for(int probe=0; probe<3; probe++){
				if(loo!=null){
					if(backgrounds==null){complete&=loo[end][probe].scoreWindow(rawEnd, values[end], probe*sites);}
					else{complete&=loo[end][probe].scoreDistantWindow(rawEnd, values[end], probe*sites);}
				}else if(backgrounds!=null){
					complete&=backgrounds[end][probe].scoreWindow(bases, rawEnd, values[end], probe*sites);
				}else{
					complete&=scoreWindow(probes[end][probe], bases, rawEnd, values[end], probe*sites);
				}
			}
			values[end][extras]=(stop-start+1)/(float)selectedConsensus.length;
			values[end][extras+1]=contigGC;
			values[end][extras+2]=identity;
			if(complete){
				validateValues(values[end]);
				usable|=bit;
			}
		}
		assert((usable&fixedMask)==0) : "Physically fixed ends cannot enter optional positional inference or training";
		return usable;
	}

	/** Uniform-table companion to the weighted window API. A missing packed word
	 * returns false; a real computed score must be finite. Writes directly into
	 * caller-owned features, preserving RrnaEndpointVector's order and arithmetic.
	 */
	private static boolean scoreWindow(final RrnaPositionalKmerTable.Table table, final byte[] bases,
			final int endpoint, final float[] output, final int from){
		assert(from>=0 && from<=output.length-table.sites) : "Each probe owns its complete resource-sized feature slice";
		for(int site=0; site<table.sites; site++){
			final long position=(long)endpoint+table.offset+site-table.radius;
			if(position<0 || position+table.k>bases.length){return false;}
			final int word=RrnaPositionalKmerTable.encode(bases, (int)position, table.k);
			if(word<0){return false;}
			final float value=(float)table.scoreKey(word);
			if(!Float.isFinite(value)){throw new IllegalArgumentException("A computed positional word score must fit the finite neural input domain");}
			output[from+site]=value;
		}
		return true;
	}
	private static void validateGC(final float value){
		if(!Float.isFinite(value) || value<0 || value>1){
			throw new IllegalArgumentException("Positional features require a measured full-contig GC fraction");
		}
	}
	private static void validateIdentity(final float value){
		if(!Float.isFinite(value) || value<0 || value>1){
			throw new IllegalArgumentException("Computed positional candidate identity must be a finite fraction");
		}
	}
	private static void validateValues(final float[] values){
		for(final float value:values){
			if(!Float.isFinite(value)){throw new IllegalArgumentException("An available positional vector cannot contain nonfinite computed features");}
		}
	}

	/** Resource-sized feature arrays,5prime then3prime; borrowed until next fill. */
	public final float[][] values;
	/** Inclusive signed correction radius, in bases. */
	public final int radius;
	private final RrnaPositionalKmerTable.Table[] tables;
	private final RrnaPositionalKmerTable.Table[][] probes;
	/** Null preserves old uniform-background train/serve arithmetic and bounds. */
	private final NcrnaWeightedEndpointTable[][] backgrounds;
	private final NcrnaPositionalLeaveOneOut[][] loo;
	private final AlignmentStats stats=new AlignmentStats(true);
	private final AlignmentStats fuzzStats;
}
