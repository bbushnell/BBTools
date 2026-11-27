package prok;

import java.util.Arrays;

import map.LongIntMap;
import prok.GradeCdsStartChanges.ReferenceContig;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.PanelAnnotationGrade.Feature;
import structures.IntList;

/** Repeated same-stop CDS start comparison against one named native control.
 * Uses GradeCdsStartChanges' source rules: a reliable changed-start pair needs
 * one selected call per stop in each arm, agreeing annotated starts, complete
 * reference flags and IDs unique across the entire original input. Overlap TP
 * is deliberately separate; a missing stop is not proof of a lost biological gene.
 * Fixed stop groups and annotation eligibility are prepared once. Trials copy
 * only selected primitive indices and reuse count/result buffers, with no sorts,
 * per-call maps or result objects. One instance belongs to one thread.
 * @author Brian, Raiden
 */
public final class NcrnaCdsStartChanges{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaCdsStartChanges(final Feature[] candidates, final ReferenceContig reference){
		require(candidates!=null && reference!=null, "Start comparison needs fixed candidates and their original source contig");
		calls=candidates.clone();
		circle=reference.circularLength();
		final int capacity=Math.addExact(calls.length, reference.rows());
		callGroups=new int[calls.length];
		truthRows=new int[capacity];
		final int[] truthStarts=new int[capacity];
		truthSpans=new int[capacity];
		truthFlags=new int[capacity];
		counts=new int[capacity];
		baselineCounts=new int[capacity];
		selectedByGroup=new int[capacity];
		seen=new boolean[calls.length];
		order=new IntList(calls.length);
		baselineOrder=new IntList(calls.length);
		final LongIntMap groups=new LongIntMap(16);
		for(int i=0; i<reference.rows(); i++){
			final Feature feature=reference.feature(i);
			final int group=group(groups, feature);
			final int start=physical(feature.locus.strand=='+' ? feature.locus.start : feature.locus.stop);
			if(truthRows[group]>0 && truthStarts[group]!=start){truthFlags[group]|=AMBIGUOUS;}
			if(!reference.completeUniqueRow(i)){truthFlags[group]|=UNCERTAIN;}
			truthRows[group]++;
			truthStarts[group]=start;
			truthSpans[group]=feature.locus.stop-feature.locus.start;
		}
		for(int i=0; i<calls.length; i++){
			final Feature feature=calls[i];
			require(feature!=null && feature.locus!=null && reference.name().equals(feature.locus.seqid)
				&& "CDS".equals(feature.locus.family) && (feature.locus.strand=='+' || feature.locus.strand=='-')
				&& feature.locus.stop<=reference.length(),
				"Native CDS candidates must belong to this exact physical reference contig and known strand");
			callGroups[i]=group(groups, feature);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Repeated trials       ----------------*/
	/*--------------------------------------------------------------*/

	/** Returns a read-only borrowed array in name(index) order. A trial invalidates
	 * previous results, including on failure. Results are zero before freezing a
	 * baseline; every later trial compares with that same control, not the last trial.
	 */
	public long[] grade(final IntList selection){
		ready=false;
		require(selection!=null, "Start comparison needs explicit native selected candidate indices");
		Arrays.fill(stats, 0);
		Arrays.fill(counts, 0);
		Arrays.fill(seen, false);
		order.clear();
		for(int i=0; i<selection.size; i++){
			final int candidate=selection.get(i);
			require(candidate>=0 && candidate<calls.length && !seen[candidate],
				"A start-comparison candidate must be distinct and inside its fixed pool");
			seen[candidate]=true;
			order.add(candidate);
			final int group=callGroups[candidate];
			counts[group]++;
			selectedByGroup[group]=candidate;
		}
		if(baseline!=null){
			for(int i=0; i<baselineOrder.size; i++){
				final int before=baselineOrder.get(i), group=callGroups[before];
				if(baselineCounts[group]>1 || counts[group]>1){stats[8]++;}
				else if(counts[group]==0){stats[1]++;}
				else{
					final Feature b=calls[before], a=calls[selectedByGroup[group]];
					if(b.locus.start!=a.locus.start || b.locus.stop!=a.locus.stop){changed(group, b, a);}
				}
			}
			for(int i=0; i<order.size; i++){
				final int group=callGroups[order.get(i)];
				if(baselineCounts[group]>1 || counts[group]>1){stats[8]++;}
				else if(baselineCounts[group]==0){stats[2]++;}
			}
			assert(stats[0]==stats[3]+stats[9]+stats[10]+stats[11])
				: "Each unique changed-start pair must be reliable, unannotated, ambiguous or uncertain";
		}
		ready=true;
		return stats;
	}
	/** Freezes one completed selected control; empty controls are valid. */
	public void freezeBaseline(final String name){
		checkResult();
		require(baseline==null && name!=null && !name.trim().isEmpty(), "Freeze one explicitly named completed CDS start control");
		System.arraycopy(counts, 0, baselineCounts, 0, counts.length);
		for(int i=0; i<order.size; i++){baselineOrder.add(order.get(i));}
		baseline=name;
	}
	public long[] statistics(){checkResult(); return stats;}
	public static String name(final int index){return NAMES[index];}
	public static int statisticsCount(){return NAMES.length;}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private void changed(final int group, final Feature before, final Feature after){
		assert(stopKey(before)==stopKey(after))
			: "A source-relative start error is a length difference only at the same strand and biological stop";
		stats[0]++;
		if(truthRows[group]==0){stats[9]++; return;}
		if((truthFlags[group]&AMBIGUOUS)!=0){stats[10]++; return;}
		if((truthFlags[group]&UNCERTAIN)!=0){stats[11]++; return;}
		stats[3]++;
		final int oldSpan=before.locus.stop-before.locus.start, newSpan=after.locus.stop-after.locus.start;
		// GradeCdsStartChanges.startDelta: the shared stop anchors the full span,
		// not the nearest half-circle, even when circumference is not a multiple of3.
		final int oldDelta=truthSpans[group]-oldSpan, newDelta=truthSpans[group]-newSpan;
		if(oldDelta==0){stats[4]++;}
		if(newDelta==0){stats[5]++;}
		if(Math.abs((long)newDelta)>Math.abs((long)oldDelta)){stats[6]++;}
		if(newSpan<oldSpan && newDelta>0){stats[7]++;}
	}
	private int group(final LongIntMap groups, final Feature feature){
		final long key=stopKey(feature);
		int group=groups.get(key);
		if(group<0){
			group=groupCount++;
			groups.put(key, group);
		}
		return group;
	}
	private long stopKey(final Feature feature){
		final Locus locus=feature.locus;
		assert(locus.strand=='+' || locus.strand=='-')
			: "Biological stop grouping requires a known transcript strand before choosing the genomic endpoint";
		return ((long)physical(locus.strand=='+' ? locus.stop : locus.start)<<1)|(locus.strand=='-' ? 1 : 0);
	}
	private int physical(final int coordinate){return circle==0 ? coordinate : (coordinate-1)%circle+1;}
	private void checkResult(){require(ready, "CDS start statistics require a completed selected trial");}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final Feature[] calls;
	private final int circle;
	private final int[] callGroups, truthRows, truthSpans, truthFlags, counts, baselineCounts, selectedByGroup;
	private final boolean[] seen;
	private final IntList order, baselineOrder;
	private final long[] stats=new long[NAMES.length];
	private int groupCount;
	private String baseline;
	private boolean ready;
	private static final int AMBIGUOUS=1, UNCERTAIN=2;
	private static final String[] NAMES={"same_stop_start_changed", "stop_not_called_rows", "added_stop_rows",
		"reliable_changed_start_pairs", "same_stop_changed_start_exact_before", "same_stop_changed_start_exact_after",
		"changed_start_error_worsened", "shortened_past_annotated_start", "ambiguous_call_stop_rows",
		"changed_start_no_same_stop_annotation", "changed_start_ambiguous_annotation", "changed_start_uncertain_annotation"};
}
