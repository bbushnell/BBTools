package prok;

import java.util.ArrayList;

import map.IntHashMap;
import map.IntHashSet;
import shared.KillSwitch;
import shared.Tools;

/** Private RNA-aware start experiment; the caller chooses which accepted families participate.
 * Suppresses an upstream start only when a retained, real downstream start at the
 * same stop/frame satisfies the existing RNA/CDS path-edge rule. Scores and RNA
 * bounds are untouched. An accepted RNA can still lose the final path; this
 * experiment therefore requires measured CDS harms, not an assumed RNA recovery.
 * @author Brian, Raiden
 */
public final class NcrnaCdsStartSelection{

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Lists contain native scored/pruned CDS candidates in genomic coordinates.
	 * Stable compaction preserves candidate identity and tie order. "Clears" permits
	 * residual overlap within the applicable cap; it does not mean zero overlap.
	 */
	public int apply(final ArrayList<Orf>[] frames, final Orf rna,
			final int sameCap, final int oppositeCap){
		if(frames==null || frames.length!=6 || rna==null || rna.type==ProkObject.CDS
				|| sameCap<0 || oppositeCap<0){
			throw new IllegalArgumentException("RNA-aware starts require six native lists, an accepted RNA and nonnegative overlap caps");
		}
		int removed=0;
		for(int f=0; f<frames.length; f++){
			final ArrayList<Orf> list=frames[f];
			if(list==null){throw new IllegalArgumentException("The retained CDS frame lists must be explicit, including empty frames");}
			final int strand=f/3;
			final int cap=strand==rna.strand ? sameCap : oppositeCap;
			blocked.clear();
			witnesses.clear();
			for(final Orf cds:list){
				if(cds.type!=ProkObject.CDS){continue;}
				assert(cds.strand==strand && cds.frame==f%3 && cds.scaflen==rna.scaflen
					&& cds.scafName.equals(rna.scafName))
					: KillSwitch.assertDie("GeneCaller frame/contig membership must survive RNA-aware start selection");
				if(blockedPrefix(cds, rna, cap)){blocked.add(stop(cds));}
			}
			if(blocked.isEmpty()){continue;}
			for(final Orf cds:list){
				if(cds.type!=ProkObject.CDS || cds.startTruncated || !blocked.contains(stop(cds))){continue;}
				final boolean legal=strand==0 ? cds.isValidPrev(rna, cap) : rna.isValidPrev(cds, cap);
				if(legal){
					final int key=stop(cds), start=start(cds), old=witnesses.get(key);
					if(old<0 || (strand==0 ? start>old : start<old)){witnesses.put(key, start);}
				}
			}
			int retained=0;
			for(int i=0, size=list.size(); i<size; i++){
				final Orf cds=list.get(i);
				boolean suppress=false;
				if(cds.type==ProkObject.CDS && blockedPrefix(cds, rna, cap)){
					final int witness=witnesses.get(stop(cds));
					final int distance=strand==0 ? witness-start(cds) : start(cds)-witness;
					suppress=witness>=0 && distance>0;
					assert(!suppress || distance%3==0)
						: KillSwitch.assertDie("Retained downstream starts at one stop/frame must differ by whole codons (GeneCaller.breakOrf)");
				}
				if(suppress){removed++;}
				else{list.set(retained++, cds);}
			}
			while(list.size()>retained){list.remove(list.size()-1);}
		}
		return removed;
	}

	/** A removable 5-prime prefix must end before the preserved biological stop.
	 * Containment alone and overlap at the 3-prime stop do not justify suppression.
	 * A contained RNA needs BOTH excess overlap and a legal retained downstream start.
	 */
	private static boolean blockedPrefix(final Orf cds, final Orf rna, final int cap){
		assert(cds.type==ProkObject.CDS && rna.type!=ProkObject.CDS)
			: "Prefix selection compares one retained coding candidate with one accepted RNA";
		final boolean beforeStop=cds.strand==0 ? rna.stop<cds.stop : rna.start>cds.start;
		return beforeStop && Tools.min(cds.stop, rna.stop)-Tools.max(cds.start, rna.start)+1>cap;
	}

	private static int stop(final Orf cds){return cds.strand==0 ? cds.stop : cds.start;}
	private static int start(final Orf cds){return cds.strand==0 ? cds.start : cds.stop;}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	// Every removed start has a strictly downstream witness. If a later RNA removes
	// that witness, its own witness is further downstream with the SAME stop/frame.
	// Moving only the 5-prime end downstream preserves every earlier legal RNA edge.
	// Thus finite stop groups retain a real witness after ALL RNA passes, in either order.
	private final IntHashSet blocked=new IntHashSet(16);
	private final IntHashMap witnesses=new IntHashMap(16);
}
