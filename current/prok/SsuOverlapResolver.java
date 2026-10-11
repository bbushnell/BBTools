package prok;

import java.util.ArrayList;
import java.util.Arrays;
import align2.PacBioScoreParameters;
import dna.AminoAcid;
import map.LongIntMap;
import map.ObjectMap;
import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;

/** One-way admission veto for clipped mitochondrial SSU rescues with a better
 * admitted 16S full-model fit. Candidate scores and boundaries are never changed.
 * Worker-local; caches are reset between contigs and strands.
 * @author Raiden */
final class SsuOverlapResolver {
	SsuOverlapResolver(boolean trace_, boolean verify_){trace=trace_;verify=verify_;}

	static boolean hasClippedMito(ArrayList<Orf>[] lists){
		assert(lists!=null):"SSU admission comparison follows RNA candidate generation";
		for(ArrayList<Orf> list:lists){if(list!=null){for(Orf o:list){if(isClippedMito(o)){return true;}}}}
		return false;
	}

	/** Lists and bases use genomic coordinates after GeneCaller.makeRnas restores the input. */
	void apply(ArrayList<Orf>[] lists, byte[] bases){
		assert(lists!=null && lists.length==2 && bases!=null):KillSwitch.assertDie("SSU arbitration requires both genomic-strand candidate lists and the restored input");
		assert(decisions.isEmpty()):KillSwitch.assertDie("Final-path audit must consume the preceding contig's SSU decisions before reuse");
		for(int strand=0; strand<2; strand++){
			for(LongIntMap map:cache.values()){if(map!=null){map.clear();}}
			final ArrayList<Orf> list=lists[strand];if(list==null){continue;}
			for(int i=0; i<list.size(); i++){
				final Orf mito=list.get(i);
				Decision winner=null;
				if(isClippedMito(mito)){
					for(Orf ssu:list){
						if(ssu==null || ssu.type!=ProkObject.r16S || ssu.ssuFullModel==null || !overlap(ssu, mito)){continue;}
						assert(ssu.strand==strand && mito.strand==strand):KillSwitch.assertDie("Genomic RNA candidates must stay in their own strand list");
						final int start=Tools.min(ssu.start, mito.start),stop=Tools.max(ssu.stop, mito.stop);
						assert(start>=0 && stop>=start && stop<bases.length):KillSwitch.assertDie("Shared SSU comparison span must fit the unchanged input contig");
						final long key=((long)start<<32)|(stop&0xffffffffL);
						byte[] span=null;
						if(!cached(ssu.ssuFullModel,key) || !cached(mito.ssuFullModel,key)){
							span=Arrays.copyOfRange(bases,start,stop+1);Tools.toUpperCase(span);
							if(strand==1){AminoAcid.reverseComplementBasesInPlace(span);}
						}
						final float ssuId=identity(ssu.ssuFullModel,key,span),mitoId=identity(mito.ssuFullModel,key,span);
						pairsCompared++;
						// Strict greater-than preserves ties. No neighboring candidate is reached transitively.
						if(ssuId>mitoId && (winner==null || ssuId>winner.ssuIdentity)){
							winner=new Decision(ssu,mito,start,stop,ssuId,mitoId);
						}
					}
				}
				if(winner!=null){list.set(i,null);decisions.add(winner);vetoCount++;}
			}
			Tools.condenseStrict(list);
		}
	}

	/** Direct reciprocal-half overlap; never forms connected components or changes a neighbor. */
	static boolean overlap(Orf a, Orf b){
		assert(a!=null && b!=null):"Overlap applies to actual admitted candidates";
		if(a.strand!=b.strand || !a.scafName.equals(b.scafName)){return false;}
		final long overlap=(long)Tools.min(a.stop,b.stop)-Tools.max(a.start,b.start)+1;
		return overlap>0 && 2*overlap>=a.length() && 2*overlap>=b.length();
	}

	private boolean cached(byte[] model, long key){
		assert(model!=null && model.length>0):"Full accepted models are required, never reconstructed from clipped alignments";
		final LongIntMap m=cache.get(model);return m!=null && m.contains(key);
	}

	private float identity(byte[] model, long key, byte[] span){
		LongIntMap m=cache.get(model);
		if(m==null){m=new LongIntMap(16);cache.put(model,m);}
		final int bits=m.get(key);if(bits>=0){return Float.intBitsToFloat(bits);}
		assert(span!=null && span.length>0):KillSwitch.assertDie("An uncached model/span fit requires the exact shared physical span");
		if(aligner==null){final PacBioScoreParameters d=PacBioScoreParameters.DEFAULT;
			aligner=new Euk18sPacBioAligner(PacBioScoreParameters.experimental(d.del4,d.del5,d.ins4,d.ins,d.del,.66f),true);}
		final float id=aligner.align(model,span,false).flatIdentity;
		assert(Float.isFinite(id) && id>=0f && id<=1f):KillSwitch.assertDie("Full-model flat identity must be finite in [0,1] before changing admission");
		alignmentCount++;alignedBases+=span.length;m.put(key,Float.floatToIntBits(id));return id;
	}

	/** Audit actual final output. Verification runs fail on neither-SSU; production does not
	 * silently add an unmeasured rollback or force a locally preferred candidate into the path. */
	void finish(ArrayList<Orf> path){
		assert(path!=null):"Final SSU outcomes require the actual shared-DP result";
		int missing=0,lost=0;
		for(Decision d:decisions){
			boolean preferred=false,anySsu=false;
			for(Orf o:path){if(o==d.ssu){preferred=true;}if(o.isSSU() && overlap(o,d.mito)){anySsu=true;}}
			if(preferred){preferredRetained++;}else{lost++;preferredLost++;}if(!anySsu){missing++;neitherRetained++;}
			if(trace || verify){
				line.clear().append("SSU_VETO\t").append(d.ssu.scafName).tab().append(d.ssu.strand)
					.tab().append(d.ssu.start+1).tab().append(d.ssu.stop+1).tab().append(d.ssu.ssuFullModelName)
					.tab().append(d.mito.start+1).tab().append(d.mito.stop+1).tab().append(d.mito.ssuFullModelName)
					.tab().append(d.start+1).tab().append(d.stop+1).tab().appendSlow(d.ssuIdentity).tab().appendSlow(d.mitoIdentity)
					.tab().append(preferred).tab().append(anySsu);
				System.err.println(line);
			}
		}
		decisions.clear();
		if(verify && missing>0){throw new IllegalStateException("Neither SSU retained after veto at "+missing+" locus/loci; see SSU_VETO rows, no rollback applied");}
		if(verify && lost>0){throw new IllegalStateException("Preferred 16S not retained after veto at "+lost+" locus/loci although another SSU survives; see SSU_VETO rows");}
	}

	private static boolean isClippedMito(Orf o){
		assert(o!=null):"RNA lists must be condensed before SSU admission comparison";
		return o.type==ProkObject.r18S && o.ssuClippedMitoRescue && o.ssuFullModel!=null
			&& Euk18sRuntimeConfig.FAMILY.equals(o.ncrnaFamily) && "euk18S_mitochondria".equals(o.ssuFullModelName);
	}

	private static final class Decision{
		Decision(Orf s, Orf m, int a, int z, float si, float mi){ssu=s;mito=m;start=a;stop=z;ssuIdentity=si;mitoIdentity=mi;}
		final Orf ssu,mito;final int start,stop;final float ssuIdentity,mitoIdentity;
	}
	final boolean trace,verify;
	long alignmentCount,alignedBases,pairsCompared,vetoCount,preferredRetained,preferredLost,neitherRetained;
	private Euk18sPacBioAligner aligner;
	private final ObjectMap<byte[],LongIntMap> cache=new ObjectMap<byte[],LongIntMap>(byte[].class,LongIntMap.class);
	private final ArrayList<Decision> decisions=new ArrayList<Decision>();
	private final ByteBuilder line=new ByteBuilder();
}
