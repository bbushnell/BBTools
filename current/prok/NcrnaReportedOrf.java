package prok;

import java.util.ArrayList;
import java.util.Collections;

import shared.KillSwitch;

/** RNA candidate whose refined endpoints are applied only after final selection.
 * Discovery, length gates, overlap resolution and path scoring use ordinary
 * start/stop fields containing the raw accepted alignment. Only glycine creates
 * this subtype; other RNA and CDS objects acquire no additional fields.
 * Stored reporting bounds are zero-based inclusive in the original sense frame.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReportedOrf extends Orf{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaReportedOrf(final String name, final int start, final int stop,
			final int strand, final byte[] bases){
		super(name, start, stop, strand, 0, bases, false, RNA);
		assert(start>=0 && stop>=start && stop<bases.length)
			: KillSwitch.assertDie("Deferred reporting requires a physical accepted RNA span");
		reportedStart=start;
		reportedStop=stop;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Saves the proposed reporting bounds and restores raw detection bounds.
	 * A proposal outside the configured length range retains the raw report;
	 * refinement cannot reject an otherwise accepted raw candidate. Physical
	 * coordinate corruption is an internal error, not a length-policy fallback.
	 */
	public void captureRefinement(final int rawStart, final int rawStop,
			final int minLen, final int maxLen){
		assert(flipped()==0 && rawStart>=0 && rawStop>=rawStart && rawStop<scaflen
			&& start>=0 && stop>=start && stop<scaflen && minLen>0 && maxLen>=minLen)
			: KillSwitch.assertDie("Capture must precede GeneCaller strand flipping and preserve physical RNA bounds");
		final long length=(long)stop-start+1;
		reportedStart=length>=minLen && length<=maxLen ? start : rawStart;
		reportedStop=length>=minLen && length<=maxLen ? stop : rawStop;
		start=rawStart;
		stop=rawStop;
	}

	/** Applies saved coordinates in the object's current frame; repeat calls are idempotent.
	 * Call only after all acceptance and selection decisions involving this object.
	 */
	public boolean applyReportedBounds(){
		assert(reportedStart>=0 && reportedStop>=reportedStart && reportedStop<scaflen)
			: KillSwitch.assertDie("Final output must use the physical bounds saved before strand flipping");
		final int nextStart=flipped()==0 ? reportedStart : scaflen-1-reportedStop;
		final int nextStop=flipped()==0 ? reportedStop : scaflen-1-reportedStart;
		final boolean changed=start!=nextStart || stop!=nextStop;
		start=nextStart;
		stop=nextStop;
		return changed;
	}

	/** Final-output handoff shared by GeneCaller and direct-scavenger consumers.
	 * Does not add/remove calls or change scores. Restore output ordering only if
	 * coordinates moved; the all-other-family path needs no sort or allocation.
	 */
	public static void applyTo(final ArrayList<Orf> calls){
		assert(calls!=null) : "Final-output handoff requires a completed candidate selection";
		boolean changed=false;
		for(Orf orf:calls){
			if(orf instanceof NcrnaReportedOrf){
				changed|=((NcrnaReportedOrf)orf).applyReportedBounds();
			}
		}
		if(changed){Collections.sort(calls);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int reportedStart, reportedStop;
}
