package prok;

import java.util.ArrayList;

import fileIO.ByteStreamWriter;
import structures.ByteBuilder;

/** Optional pre-report selected-path evidence. Borrows the caller's native writer.
 * All coordinates are genomic, zero-based inclusive. P rows preserve every selected
 * feature, including multiplicity and path bounds; their local indices are not
 * across-run identities. Q rows compare the same objects in both geometries.
 * No scoring, selection, source mutation or sequence retention occurs here.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSelectionTrace{

	/** Called before report application, while predecessor links still describe the selected chain.
	 * Target-vs-all inspection costs O(R*N), for R moved target RNAs and N selected
	 * features. It creates no per-pair objects and emits each unordered pair once.
	 */
	public static void write(final ByteStreamWriter output, final ByteBuilder b, final int ordinal,
			final long readID, final String contig, final String target, final ArrayList<Orf> path,
			final int sameCap, final int oppositeCap){
		if(output==null){return;}
		require(b!=null && ordinal>=0 && readID>=0 && path!=null && sameCap>=0 && oppositeCap>=0,
			"Selection evidence requires the actual input/contig, completed path and caller overlap limits");
		b.clear();
		int pairs=0;
		for(int i=0; i<path.size(); i++){
			final Orf orf=path.get(i);
			require(orf.start>=0 && orf.stop>=orf.start && orf.stop<orf.scaflen,
				"Selected path evidence must precede physical report application");
			prefix(b, 'P', ordinal, readID, contig).tab().append(i).tab().append(orf.strand).tab()
				.append(orf.type).tab().append(orf.frame).tab().append(orf.ncrnaFamily==null ? "." : orf.ncrnaFamily).tab()
				.append(orf.trnaModel==null ? "." : orf.trnaModel).tab()
				.append(orf.start).tab().append(orf.stop).tab().append(reportStart(orf)).tab().append(reportStop(orf)).tab()
				.append(Float.floatToIntBits(orf.orfScore)).nl();
		}
		for(int i=0; i<path.size(); i++){
			final Orf a=path.get(i);
			if(!movedTarget(a, target)){continue;}
			for(int j=0; j<path.size(); j++){
				final Orf c=path.get(j);
				if(i==j || j<i && movedTarget(c, target)){continue;}
				final boolean edge=a.prev()==c || c.prev()==a;
				final int before=overlap(a.start, a.stop, c.start, c.stop);
				final int after=overlap(reportStart(a), reportStop(a), reportStart(c), reportStop(c));
				if(!edge && before==0 && after==0){continue;}
				pairs++;
				prefix(b, 'Q', ordinal, readID, contig).tab().append(i).tab().append(j).tab().append(edge ? 1 : 0).tab()
					.append(a.strand==c.strand ? 1 : 0).tab().append(c.type).tab().append(before).tab().append(after).tab()
					.append(a.strand==c.strand ? sameCap : oppositeCap).tab()
					.append(order(a.start, a.stop, c.start, c.stop)).tab()
					.append(order(reportStart(a), reportStop(a), reportStart(c), reportStop(c))).tab()
					.append(containment(a.start, a.stop, c.start, c.stop)).tab()
					.append(containment(reportStart(a), reportStop(a), reportStart(c), reportStop(c))).nl();
			}
		}
		prefix(b, 'T', ordinal, readID, contig).tab().append(path.size()).tab().append(pairs).nl();
		synchronized(output){output.print(b);}
	}

	private static ByteBuilder prefix(final ByteBuilder b, final char kind, final int ordinal,
			final long readID, final String contig){return b.append(kind).tab().append(ordinal).tab().append(readID).tab().append(contig);}
	private static int reportStart(final Orf orf){return orf instanceof NcrnaReportedOrf ? ((NcrnaReportedOrf)orf).reportStart() : orf.start;}
	private static int reportStop(final Orf orf){return orf instanceof NcrnaReportedOrf ? ((NcrnaReportedOrf)orf).reportStop() : orf.stop;}
	private static boolean movedTarget(final Orf orf, final String target){
		return orf.ncrnaFamily!=null && (orf.ncrnaFamily.equals(target) || orf.ncrnaFamily.length()>target.length()
			&& orf.ncrnaFamily.startsWith(target) && orf.ncrnaFamily.charAt(target.length())=='_')
			&& (orf.start!=reportStart(orf) || orf.stop!=reportStop(orf));
	}
	private static int overlap(final int a, final int b, final int c, final int d){return Math.max(0, Math.min(b, d)-Math.max(a, c)+1);}
	/** Same stop-then-start order used by PFeature, including exact ties. */
	private static int order(final int a, final int b, final int c, final int d){return b==d ? Integer.compare(a, c) : Integer.compare(b, d);}
	/** 0=neither contains;1=first contains;2=second contains;3=equal spans. */
	private static int containment(final int a, final int b, final int c, final int d){return (a<=c && b>=d ? 1 : 0)|(c<=a && d>=b ? 2 : 0);}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalStateException(why);}}

	public static final String HEADER="#ncrna_selected_geometry\t1\n"
		+"#Coordinates: zero-based inclusive genomic. Compare path identity with multiplicity; report endpoints intentionally differ.\n"
		+"#I\tinput_ordinal0\tassembly\n"
		+"#P\tinput_ordinal0\tread_id0\tcontig\tlocal_index0\tstrand\ttype\tframe\tfamily\tmodel\tpath_start0\tpath_stop0\treport_start0\treport_stop0\tscore_bits\n"
		+"#Q\tinput_ordinal0\tread_id0\tcontig\tfirst_index0\tsecond_index0\tselected_edge\tsame_strand\tsecond_type\tpath_overlap\treport_overlap\toverlap_cap\tpath_order\treport_order\tpath_containment\treport_containment\n"
		+"#T\tinput_ordinal0\tread_id0\tcontig\tselected_features\tpair_rows\n";
}
