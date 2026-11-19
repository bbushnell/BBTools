package benchmark;

import java.util.ArrayList;
import java.util.List;

import gff.GffLine;
import prok.Orf;
import prok.ProkObject;
import structures.ByteBuilder;
import structures.IntList;

/** Thread-confined compact calls retained across the calling/grading barrier.
 * Four primitive integers per CDS; no caller, sequence or classifier evidence
 * survives. Only fields consumed by the accepted grader are reconstructed.
 * The owner appends serially, then reads after caller workers have joined.
 * Reconstructed rows are grading views, not complete copies of emitted GFF.
 * @author Sayu, Fischl
 */
public final class CallGenesLeanCalls{

	/*--------------------------------------------------------------*/
	/*----------------     Capture and grading      ----------------*/
	/*--------------------------------------------------------------*/

	/** Append a borrowed contig batch in input order, preserving duplicate CDS and partial bits.
	 * Null or empty batches add nothing. Every feature contributes to features(),
	 * but only CDS coordinates are retained, converted to one-based inclusive GFF bounds.
	 * The contig key must match Orf.appendGff name truncation. Failed appends may
	 * leave partial state; discard this capture after invalid input. */
	public void append(String contig, List<Orf> calls){
		require(contig!=null && !contig.isEmpty() && contig.indexOf(' ')<0 && contig.indexOf('\t')<0,
			"Compact calls require the exact emitted GFF contig key");
		if(calls==null || calls.isEmpty()){return;}
		final int id=contigs.size();
		contigs.add(contig);
		features=Math.addExact(features, calls.size());
		for(Orf call:calls){
			require(call!=null, "A missing called feature cannot be omitted from conservation");
			if(call.type!=ProkObject.CDS){continue;}
			require(matchesContig(call.scafName, contig) && call.start>=0 && call.stop>=call.start
				&& (call.strand==0 || call.strand==1), "Invalid compact CDS coordinates or contig identity");
			data.add(id);
			data.add(Math.addExact(call.start, 1));
			data.add(Math.addExact(call.stop, 1));
			// Orf.appendGff uses these exact forward-coordinate partial bits and phase0.
			data.add(call.strand|(call.partialMask()<<1));
		}
	}

	/** Reconstruct one genome/arm for grading after all caller workers have joined.
	 * Returns fresh records and a fresh list; the retained capture is unchanged.
	 * GffLine(Orf) omits emitted partial metadata, so it cannot establish grade parity. */
	public ArrayList<GffLine> lines(){
		require(data.size%4==0 && GffLine.parseSeqid && GffLine.parseType && GffLine.parseAttributes,
			"Compact grading requires intact records and the accepted GFF parser fields");
		final ArrayList<GffLine> result=new ArrayList<GffLine>(size());
		final ByteBuilder row=new ByteBuilder();
		for(int i=0; i<data.size; i+=4){
			final int packed=data.array[i+3], partial=packed>>>1;
			row.clear().append(contigs.get(data.array[i])).append("\tBBTools\tCDS\t")
				.append(data.array[i+1]).tab().append(data.array[i+2]).append("\t.\t")
				.append((packed&1)==0 ? '+' : '-').append("\t0\t");
			if(partial==0){
				row.append('.');
			}else{
				row.append("partial=").append((partial&1)==0 ? '0' : '1').append((partial&2)==0 ? '0' : '1');
			}
			result.add(new GffLine(row.toBytes()));
		}
		return result;
	}
	/** Number of retained CDS rows, including duplicates. */
	public int size(){return data.size/4;}
	/** Total appended features before the CDS filter, including RNA. */
	public long features(){return features;}

	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/

	/** Match the name truncation performed by Orf.appendGff without per-gene substrings. */
	private static boolean matchesContig(String name, String contig){
		return name!=null && name.startsWith(contig) && (name.length()==contig.length()
			|| name.charAt(contig.length())==' ' || name.charAt(contig.length())=='\t');
	}
	/** Reject invalid capture state even when assertions are disabled. */
	private static void require(boolean value, String message){
		if(!value){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Owned flat records: contig ordinal, GFF start, GFF stop, and packed strand/partial bits. */
	private final IntList data=new IntList();
	/** One emitted contig key per nonempty appended batch; record ordinals index this list. */
	private final ArrayList<String> contigs=new ArrayList<String>();
	/** Count before CDS filtering; updated only by the capture owner. */
	private long features;
}
