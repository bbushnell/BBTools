package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Objects;

import map.ObjectIntMap;
import shared.Tools;
import structures.ByteBuilder;

/** One complete retained candidate graph and its same-run path oracle.
 * Observation copies scalar metadata into isolated nodes, never native objects
 * or DNA. Repeated trials reuse those nodes and the original six list orders.
 * Caller policy/PGM binding and saved-file transport belong to the driver.
 * @author Brian, Raiden
 */
public final class NcrnaPathSnapshot implements GeneCallerPathTrace.Observer{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaPathSnapshot(final String source_, final long read_, final String contig_){
		require(source_!=null && !source_.isEmpty() && read_>=0 && contig_!=null && !contig_.isEmpty(),
			"A graph needs explicit source/read/contig identity, including empty graphs");
		source=source_;
		read=read_;
		contig=contig_;
	}

	@Override
	@SuppressWarnings("unchecked")
	public void begin(final String name, final int length, final int count){
		require(nodes==null && length>=0 && count>=0 && (count==0 ? name==null : contig.equals(name)),
			"Capture exactly one complete graph for the declared contig");
		contigLength=length;
		nodes=new Orf[count];
		ids=new long[count];
		initial=new int[Math.multiplyExact(count, WIDTH)];
		expected=new int[Math.multiplyExact(count, 4)];
		links=new long[Math.multiplyExact(count, 2)];
		labels=new String[Math.multiplyExact(count, 3)];
		seen=new boolean[count];
		frames=new ArrayList[6];
		for(int f=0; f<6; f++){frames[f]=new ArrayList<Orf>();}
	}

	@Override
	public boolean candidate(final long id, final int start, final int stop, final int strand,
			final int frame, final int type, final int stopLeft, final float score,
			final float bonus, final boolean nativeMeasured, final float nativeScore){
		require(nodes!=null && !ended && !pending && received<nodes.length,
			"Each retained candidate needs one immediately following metadata record");
		require(strand>=0 && strand<2 && frame>=0 && frame<3 && !nativeMeasured && bonus==0,
			"This graph format excludes protein/cascade scoring dependencies");
		final int list=3*strand+frame;
		require(id==pack(list, frames[list].size()) && list>=lastList && Float.isFinite(score),
			"Capture dense original frame-list identities and finite scores in native list order");
		pending=true;
		pendingId=id;
		pendingStart=start;
		pendingStop=stop;
		pendingType=type;
		pendingScore=score;
		lastList=list;
		return false;// Candidate and state records suffice; do not produce edge/window rows.
	}

	@Override
	public void nodeClass(final long id, final String className){
		require(pending && pendingId==id && (Orf.class.getName().equals(className)
			|| NcrnaReportedOrf.class.getName().equals(className)),
			"Path replay does not support unrecorded candidate-subtype state: "+className);
	}

	@Override
	public void metadata(final long id, final int flipped, final int startCodon, final int stopCodon,
			final boolean startTruncated, final boolean stopTruncated, final float startScore,
			final float stopScore, final float kmerScore, final float oppositePenalty,
			final String model, final String anticodon, final String family, final int partial,
			final float acceptedRrnaIdentity, final boolean deferred, final int reportedStart,
			final int reportedStop, final boolean proteinEvidence, final boolean ssuEvidence){
		require(pending && pendingId==id && !proteinEvidence && !ssuEvidence && oppositePenalty==0,
			"Replay requires the matching node metadata without unsupported protein/SSU/opposite-penalty dependencies");
		require(Float.isFinite(startScore) && Float.isFinite(stopScore) && Float.isFinite(kmerScore)
			&& Float.isFinite(acceptedRrnaIdentity) && partial>=0 && partial<=3,
			"Captured report and overlap-rescoring components must remain finite and physically meaningful");
		final int list=(int)(id>>>32);
		final int strand=list/3;
		final int frame=list%3;
		require(!deferred || pendingType==ProkObject.RNA && frame==0,
			"Deferred report nodes retain their native RNA type and placeholder frame");
		require(deferred || reportedStart==-1 && reportedStop==-1,
			"Only deferred nodes may carry independent report bounds");
		final Orf node=deferred ? new NcrnaReportedOrf(contig, pendingStart, pendingStop, strand,
			contigLength, flipped, startCodon, stopCodon, startTruncated, stopTruncated, reportedStart, reportedStop)
			: new Orf(contig, pendingStart, pendingStop, strand, frame, contigLength, flipped, pendingType,
				startCodon, stopCodon, startTruncated, stopTruncated);
		node.startScore=startScore;
		node.stopScore=stopScore;
		node.kmerScore=kmerScore;
		node.orfScore=pendingScore;
		node.trnaModel=model;
		node.trnaAnticodon=anticodon;
		node.ncrnaFamily=family;
		node.ncrnaPartial=partial;
		node.acceptedRrnaIdentity=acceptedRrnaIdentity;
		if(!frames[list].isEmpty()){
			require(frames[list].get(frames[list].size()-1).compareTo(node)<=0,
				"Native DP list advancement depends on stable endpoint order, including exact-coordinate ties");
		}
		frames[list].add(node);
		nodes[received]=node;
		ids[received]=id;
		indices.put(node, received);
		final int at=received*WIDTH;
		initial[at+START]=node.start;
		initial[at+STOP]=node.stop;
		initial[at+SCORE]=bits(node.orfScore);
		initial[at+START_SCORE]=bits(startScore);
		initial[at+STOP_SCORE]=bits(stopScore);
		initial[at+KMER_SCORE]=bits(kmerScore);
		initial[at+FLIPPED]=flipped;
		initial[at+PARTIAL]=partial;
		initial[at+RRNA_ID]=bits(acceptedRrnaIdentity);
		labels[3*received]=model;
		labels[3*received+1]=anticodon;
		labels[3*received+2]=family;
		received++;
		pending=false;
	}

	@Override
	public void scored(final long id, final float plusScore, final float minusScore,
			final int plusLength, final int minusLength, final long plusPrev, final long minusPrev){
		require(!ended && !pending && received==nodes.length, "Scoring begins only after the complete candidate graph");
		final int index=index(id);
		require(!seen[index] && Float.isFinite(plusScore) && Float.isFinite(minusScore)
			&& plusLength>0 && minusLength>0 && plusPrev>=-1 && minusPrev>=-1,
			"Record each native predecessor state exactly once with -1 as the only absent-node sentinel");
		if(plusPrev>=0){index(plusPrev);}
		if(minusPrev>=0){index(minusPrev);}
		expected[4*index]=bits(plusScore);
		expected[4*index+1]=bits(minusScore);
		expected[4*index+2]=plusLength;
		expected[4*index+3]=minusLength;
		links[2*index]=plusPrev;
		links[2*index+1]=minusPrev;
		seen[index]=true;
		scored++;
	}

	@Override
	public void end(final long bestNode, final float bestScore, final int selectedNodes){
		require(!ended && !pending && received==nodes.length && scored==nodes.length,
			"A complete graph must retain every candidate and both scored predecessor states");
		require(selectedNodes>=0 && selectedNodes<=nodes.length && Float.isFinite(bestScore)
			&& (nodes.length==0 ? bestNode==-1 && selectedNodes==0 && bestScore==0 : bestNode>=0),
			"The native endpoint and selected count must describe the same complete graph");
		best=bestNode;
		if(best>=0){
			final int index=index(best);
			require(bits(bestScore)==bits(Tools.max(value(expected[4*index]), value(expected[4*index+1]))),
				"The native endpoint score must equal its recorded predecessor-state maximum");
		}
		selected=selectedNodes;
		ended=true;
	}

	/** Validate complete reported rows against the actual same-run native output.
	 * Invoke after the caller returns, never from an observation callback.
	 */
	public void acceptBaseline(final ArrayList<Orf> nativePath){
		require(ended && !accepted && nativePath!=null && nativePath.size()==selected,
			"Bind one completed native report before any replay trial");
		appendRows(nativePath, baselineRows);
		acceptBaselineRowsInner();
	}

	/** Saved reports retain exact bytes, not reconstructed rounded-score Orfs. */
	public void acceptBaselineRows(final ByteBuilder rows, final int count){
		require(ended && !accepted && rows!=null && count==selected,
			"A saved native report must conserve the original selected count");
		baselineRows.clear();
		baselineRows.append(rows);
		acceptBaselineRowsInner();
	}
	private void acceptBaselineRowsInner(){
		final ArrayList<Orf> copiedPath=new ArrayList<Orf>(selected);
		long next=best;
		while(next>=0){
			require(copiedPath.size()<selected, "Native predecessor chain cannot cycle or exceed its selected count");
			final int index=index(next);
			copiedPath.add(nodes[index]);
			final boolean plus=value(expected[4*index])>=value(expected[4*index+1]);
			next=links[2*index+(plus ? 0 : 1)];
		}
		require(copiedPath.size()==selected, "The saved predecessor chain must conserve selected candidates");
		Collections.sort(copiedPath);
		NcrnaReportedOrf.applyTo(copiedPath);
		selectedIds=new long[selected];
		for(int i=0; i<selected; i++){selectedIds[i]=nodeId(copiedPath.get(i));}
		appendRows(copiedPath, trialRows);
		require(sameRows(), "Copied metadata must reproduce every native feature row and attribute before replay");
		accepted=true;
	}

	/*--------------------------------------------------------------*/
	/*----------------        Repeated trials       ----------------*/
	/*--------------------------------------------------------------*/

	/** Reset all mutable scoring/report state on the same isolated candidate objects. */
	public void restore(){
		require(accepted, "Replay begins only after complete same-run baseline report validation");
		for(int i=0; i<nodes.length; i++){
			final Orf node=nodes[i];
			final int at=i*WIDTH;
			if(node.flipped()!=initial[at+FLIPPED]){node.flip();}
			node.start=initial[at+START];
			node.stop=initial[at+STOP];
			node.orfScore=value(initial[at+SCORE]);
			node.startScore=value(initial[at+START_SCORE]);
			node.stopScore=value(initial[at+STOP_SCORE]);
			node.kmerScore=value(initial[at+KMER_SCORE]);
			node.ncrnaPartial=initial[at+PARTIAL];
			node.acceptedRrnaIdentity=value(initial[at+RRNA_ID]);
			node.trnaModel=labels[3*i];
			node.trnaAnticodon=labels[3*i+1];
			node.ncrnaFamily=labels[3*i+2];
			node.proteinBonus=0;
			node.tailCascadeOppositePenalty=0;
			node.prevPlus=null;
			node.prevMinus=null;
			node.pathScorePlus=0;
			node.pathScoreMinus=0;
			node.pathLengthPlus=1;
			node.pathLengthMinus=1;
		}
		ready=true;
	}
	public int candidates(){return nodes==null ? 0 : nodes.length;}
	/** Physical captured length, including an accepted graph with no candidates. */
	public int length(){
		require(accepted, "Captured length lookup requires a complete accepted snapshot");
		return contigLength;
	}

	/** Stable frame-list identity of an isolated replay node; null is the absent
	 * predecessor sentinel -1. Native or foreign snapshot objects are rejected.
	 * Reporting sorts and deferred coordinate changes do not change this identity.
	 */
	public long candidateId(final Orf node){
		require(accepted, "Candidate joins require a complete accepted snapshot");
		return nodeId(node);
	}

	/** Exact captured local score, independent of the current trial or reporting
	 * bounds. The saved ID must belong to this complete accepted snapshot.
	 */
	public float originalScore(final long id){
		require(accepted, "Original score joins require a complete accepted snapshot");
		return value(initial[index(id)*WIDTH+SCORE]);
	}

	/** Receives potential report geometry for EVERY retained candidate, including
	 * unselected nodes, in original frame-list order. Bounds are zero-based,
	 * inclusive and forward-genomic; model/family are the nullable captured labels.
	 * Source/read/contig identity belongs to the enclosing snapshot. The receiver
	 * must not mutate this snapshot during visitation. No mutable Orf is exposed.
	 */
	public interface ReportVisitor{
		void candidate(long id, int type, int strand, int start, int stop, String model, String family);
	}

	/** Visits immutable original reporting metadata without restoring or consuming
	 * a trial. Ordinary captured bounds are already genomic, regardless of flip
	 * parity. Deferred RNA alone needs its original-sense reporting transform.
	 * No result object is allocated per candidate; callers may build a grade pool
	 * once and reuse it across score trials. This is not a selected-call report.
	 */
	public void visitCandidateReports(final ReportVisitor visitor){
		require(accepted && visitor!=null, "Report projection requires a complete accepted snapshot and receiver");
		for(int i=0; i<nodes.length; i++){visitReport(i, visitor);}
	}

	/** Visits only the accepted original selected calls, in their native report
	 * order. Uses the saved predecessor-chain identities and frozen report metadata,
	 * not current trial state. Does not restore, consume or mutate a score trial.
	 * The receiver must not mutate this snapshot. No mutable Orf is exposed.
	 */
	public void visitBaselineReports(final ReportVisitor visitor){
		require(accepted && visitor!=null, "Baseline report projection requires a complete accepted snapshot and receiver");
		for(final long id:selectedIds){visitReport(index(id), visitor);}
	}
	/** Single genomic projection shared by the full pool and saved selection. */
	private void visitReport(final int i, final ReportVisitor visitor){
		assert(accepted && i>=0 && i<nodes.length && visitor!=null)
			: "Accepted report visitors resolve only retained candidate indices";
		final Orf node=nodes[i];
		final int at=i*WIDTH;
		int start=initial[at+START], stop=initial[at+STOP];
		if(node instanceof NcrnaReportedOrf){
			final NcrnaReportedOrf deferred=(NcrnaReportedOrf)node;
			final int senseStart=deferred.reportedSenseStart(), senseStop=deferred.reportedSenseStop();
			final boolean flipped=initial[at+FLIPPED]!=0;
			start=flipped ? contigLength-1-senseStop : senseStart;
			stop=flipped ? contigLength-1-senseStart : senseStop;
		}
		assert(start>=0 && stop>=start && stop<contigLength)
			: "Saved physical Orf/report bounds must stay physical after projection: source="+source+
				" read="+read+" id="+ids[i]+" bounds="+start+".."+stop+" length="+contigLength;
		visitor.candidate(ids[i], node.type, node.strand, start, stop, labels[3*i], labels[3*i+2]);
	}

	/** Changes only an RNA contribution after restoration, never candidate geometry. */
	public void setScore(final long id, final float score){
		require(ready && Float.isFinite(score), "Trial scores must be finite and assigned after restoration");
		final Orf node=nodes[index(id)];
		require(node.type!=ProkObject.CDS, "This score experiment varies RNA while retaining native CDS sub-scores");
		node.orfScore=score;
	}

	/** Execute the existing native path search on a separately configured replay caller. */
	public ArrayList<Orf> select(final GeneCaller caller){
		require(ready && caller!=null, "Each path trial requires explicit restoration and an independent replay caller");
		ready=false;
		final ArrayList<Orf> result=caller.findPath(frames, contigLength);
		NcrnaReportedOrf.applyTo(result);
		return result;
	}

	/** Exact zero-change oracle: both states of every node, final IDs and complete rows. */
	public void verifyBaseline(final GeneCaller caller){
		restore();
		final ArrayList<Orf> result=select(caller);
		verifyAgainst(this, result);
	}

	/** Compare a completed trial to a separately observed native execution.
	 * Finite changed-score fixtures use this oracle without copying native DP code.
	 */
	public void verifyAgainst(final NcrnaPathSnapshot reference, final ArrayList<Orf> result){
		require(accepted && reference!=null && reference.accepted && result!=null
			&& nodes.length==reference.nodes.length && source.equals(reference.source)
			&& read==reference.read && contig.equals(reference.contig) && contigLength==reference.contigLength,
			"Path-state parity requires the same complete source-qualified candidate graph");
		for(int i=0; i<nodes.length; i++){
			final Orf node=nodes[i];
			require(sameCandidate(i, reference) && bits(node.orfScore)==reference.initial[i*WIDTH+SCORE],
				"Native trial and replay must retain the same candidates and exact assigned local scores");
			require(bits(node.pathScorePlus)==reference.expected[4*i]
				&& bits(node.pathScoreMinus)==reference.expected[4*i+1]
				&& node.pathLengthPlus==reference.expected[4*i+2] && node.pathLengthMinus==reference.expected[4*i+3]
				&& nodeId(node.prevPlus)==reference.links[2*i] && nodeId(node.prevMinus)==reference.links[2*i+1],
				"Replay changed a native node score bit, predecessor or path length at "+source+"/"+read+"/"+ids[i]);
		}
		require(result.size()==reference.selectedIds.length, "Replay must retain the reference execution's selected count");
		for(int i=0; i<result.size(); i++){
			require(nodeId(result.get(i))==reference.selectedIds[i], "Replay must preserve selected identity and stable final order");
		}
		appendRows(result, trialRows);
		require(sameRows(reference.baselineRows, trialRows), "Replay must preserve complete emitted rows, including deferred bounds and attributes");
	}

	private boolean sameCandidate(final int i, final NcrnaPathSnapshot other){
		final Orf a=nodes[i];
		final Orf b=other.nodes[i];
		if(ids[i]!=other.ids[i] || a.type!=b.type || a.frame!=b.frame || a.strand!=b.strand
			|| a.startCodon!=b.startCodon || a.stopCodon!=b.stopCodon
			|| a.startTruncated!=b.startTruncated || a.stopTruncated!=b.stopTruncated){return false;}
		for(int field=0; field<WIDTH; field++){
			if(field!=SCORE && initial[i*WIDTH+field]!=other.initial[i*WIDTH+field]){return false;}
		}
		for(int field=0; field<3; field++){
			if(!Objects.equals(labels[3*i+field], other.labels[3*i+field])){return false;}
		}
		if(a instanceof NcrnaReportedOrf || b instanceof NcrnaReportedOrf){
			if(!(a instanceof NcrnaReportedOrf && b instanceof NcrnaReportedOrf)){return false;}
			final NcrnaReportedOrf x=(NcrnaReportedOrf)a;
			final NcrnaReportedOrf y=(NcrnaReportedOrf)b;
			return x.reportedSenseStart()==y.reportedSenseStart() && x.reportedSenseStop()==y.reportedSenseStop();
		}
		return true;
	}

	private int index(final long id){
		final int list=(int)(id>>>32);
		final int item=(int)id;
		require(id>=0 && list>=0 && list<6 && item>=0 && item<frames[list].size(), "A saved predecessor must identify a retained frame-list node");
		return indices.get(frames[list].get(item));
	}
	private long nodeId(final Orf node){
		if(node==null){return -1;}
		final int index=indices.get(node);
		require(index>=0, "Replay predecessor objects cannot escape the isolated candidate graph");
		return ids[index];
	}
	private static void appendRows(final ArrayList<Orf> path, final ByteBuilder rows){
		assert(path!=null && rows!=null) : "Report parity compares two completed feature lists";
		rows.clear();
		for(final Orf node:path){node.appendGff(rows).nl();}
	}
	private boolean sameRows(){
		return sameRows(baselineRows, trialRows);
	}
	private static boolean sameRows(final ByteBuilder a, final ByteBuilder b){
		if(a.length()!=b.length()){return false;}
		for(int i=0; i<a.length(); i++){
			if(a.array[i]!=b.array[i]){return false;}
		}
		return true;
	}
	private static long pack(final int list, final int index){return ((long)list<<32)|(index&0xffffffffL);}
	private static int bits(final float value){return Float.floatToRawIntBits(value);}
	private static float value(final int bits){return Float.intBitsToFloat(bits);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	@Override
	public void edge(long node, long predecessor, int predecessorStrand, int overlap, float predecessorScore,
			int predecessorLength, float overlapScore, float beforeTransitions, float proposedScore, int proposedLength, boolean replaces){}
	@Override
	public void window(long node, int predecessorList, int lastPosition, int bestPosition, int nominalMinimum, int lastVisited){}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final String source, contig;
	public final long read;
	private final ObjectIntMap<Orf> indices=new ObjectIntMap<Orf>(Orf.class);
	private final ByteBuilder baselineRows=new ByteBuilder(), trialRows=new ByteBuilder();
	private ArrayList<Orf>[] frames;
	private Orf[] nodes;
	private long[] ids, links, selectedIds;
	private int[] initial, expected;
	private String[] labels;
	private boolean[] seen;
	private int contigLength, received, scored, lastList, selected;
	private boolean pending, ended, accepted, ready;
	private long pendingId, best;
	private int pendingStart, pendingStop, pendingType;
	private float pendingScore;
	private static final int START=0, STOP=1, SCORE=2, START_SCORE=3, STOP_SCORE=4,
		KMER_SCORE=5, FLIPPED=6, PARTIAL=7, RRNA_ID=8, WIDTH=9;
}
