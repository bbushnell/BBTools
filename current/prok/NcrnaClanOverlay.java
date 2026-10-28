package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Materializes a clan evidence overlay and preserves unresolved contexts for followup.
 * Own full support retains the source subtype; other-only support is clan-only.
 * Length bands are diagnostic, never a veto on independently supported cores.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaClanOverlay {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaClanOverlay.class, false);
		String joined=null, contexts=null, prefix=null, family=null, rescue=null, phase="FIRST_PASS";
		int expected=-1, low=-1, high=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use explicit flag=value arguments");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("joined")){joined=value;}else if(key.equals("contexts")){contexts=value;}
			else if(key.equals("outprefix")){prefix=value;}else if(key.equals("sourcefamily")){family=value;}
			else if(key.equals("rescue")){rescue=value;}else if(key.equals("evidencestage")){phase=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}
			else if(key.equals("minlen")){low=Integer.parseInt(value);}else if(key.equals("maxlen")){high=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(joined!=null && contexts!=null && prefix!=null && family!=null && expected>0 && low>=17 && high>=low,
			"Bind the complete joined evidence, original contexts, source family, denominator and diagnostic length band");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		require(phase.equals("FIRST_PASS") || phase.equals("RAW_ADJUDICATED"), "Declare the evidence stage without overwriting source provenance");
		try{
			final NcrnaClanOverlay x=new NcrnaClanOverlay(family, low, high, phase);
			x.load(joined, expected, false);if(rescue!=null){x.load(rescue, x.unresolved, true);}x.write(contexts, prefix);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaClanOverlay(String family_, int low_, int high_, String phase_){
		family=family_;low=low_;high=high_;phase=phase_;int selected=-1;
		for(int i=0; i<MODELS.length; i++){if(MODELS[i].equals(family)){selected=i;}}
		require(selected>=0, "The overlay requires a declared B12 clan source family");own=selected;
	}
	private void load(String path, int expected, boolean rescue){
		assert(expected>=0):"A rescue join must cover exactly the original unresolved population";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(in, header(), path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==72 && p.termEquals(family, 21), "Source family or complete three-model join schema changed");
				final String id=p.parseString(0);require(id.startsWith(family+"_"), "Source ID does not belong to its declared family: "+id);
				final byte[] source=prefix(row, 21);Entry e=byId.get(id);
				if(!rescue){require(e==null, "Duplicate original source ID: "+id);e=new Entry(id, source, p);byId.put(id, e);entries.add(e);}
				else{require(e!=null && e.initialUnresolved && !e.rescueSeen && Arrays.equals(e.source, source), "Rescue must preserve each originally unresolved source exactly once: "+id);e.rescueSeen=true;}
				classify(e, p);e.stage=rescue?"SENSITIVE_FOLLOWUP":phase;
				if(!rescue){e.initialUnresolved=e.model<0;if(e.initialUnresolved){unresolved++;}}
				count++;
			}
		}finally{require(!in.close(), "Evidence join read failed: "+path);}
		require(count==expected, "Evidence must conserve the complete declared population: "+count+" expected="+expected);
		if(rescue){for(Entry e:entries){require(!e.initialUnresolved || e.rescueSeen, "Sensitive followup cannot silently omit a requested locus: "+e.id);}}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void classify(Entry e, LineParser1 p){
		int count=0, sole=-1, coordinate=-1, bestOverlap=-1;boolean tied=false;
		final ByteBuilder support=new ByteBuilder();boolean ownFull=false;
		for(int m=0; m<MODELS.length; m++){
			final int c=26+15*m;
			final boolean any=p.termEquals("CM_FULL_SPAN_ANY_HIT_SUPPORTED", c);
			if(!any && !p.termEquals("CM_FULL_SPAN_SUPPORTED", c)){continue;}
			require((any?p.parseInt(c+2)>=1:p.parseInt(c+2)==1) && p.parseInt(c+3)==1 && p.parseInt(c+4)==CLEN[m]
				&& p.termEquals("no", c+7) && p.parseDouble(c+8)>=GA[m], "Admitted support must retain the pinned Rfam15.1 full-span/GA criteria: "+e.id);
			require(p.termEquals("EXACT", c+12) || p.termEquals("HIT_WITHIN_ANNOTATION", c+12) || p.termEquals("ANNOTATION_WITHIN_HIT", c+12)
				|| (any && p.termEquals("PARTIAL_OVERLAP", c+12)), "Unexpected source-hit geometry: "+e.id);
			final int overlap=Math.min(e.left+e.length, p.parseInt(c+6))-Math.max(e.left+1, p.parseInt(c+5))+1;
			require(overlap>0 && overlap==p.parseInt(c+9), "A nonoverlapping neighbor cannot establish source-locus admission: "+e.id);
			if(overlap>bestOverlap){coordinate=m;bestOverlap=overlap;tied=false;}
			else if(overlap==bestOverlap && (p.parseInt(c+5)!=p.parseInt(26+15*coordinate+5) || p.parseInt(c+6)!=p.parseInt(26+15*coordinate+6))){tied=true;}
			if(count>0){support.append(',');}support.append(MODELS[m]);count++;sole=m;ownFull|=m==own;
		}
		require(count==p.parseInt(22) && p.termEquals(ownFull?"true":"false", 24)
			&& p.termEquals(count==1?MODELS[sole]:"NA", 23), "Joined support fields disagree with their raw model columns: "+e.id);
		// Own support chooses the own model. Cross-labeled loci stay clan-only even
		// when maximum source overlap selects coordinates from one of several CMs.
		require(ownFull || !tied, "Equal-overlap other-only full cores require explicit coordinate adjudication: "+e.id);
		e.model=ownFull?own:coordinate;e.support=count==0?"NONE":support.toString();
		e.label=ownFull?family:count>0?"CROSS_LABELED":"UNRESOLVED";
		if(e.model<0){return;}
		final int c=26+15*e.model;e.from=p.parseInt(c+5);e.to=p.parseInt(c+6);
		e.start=p.parseInt(c+13);e.stop=p.parseInt(c+14);
		require(e.from>0 && e.to>=e.from && e.to<=e.length+e.left+e.right && e.start>0 && e.stop<=e.contigLength
			&& e.stop-e.start==e.to-e.from, "Supported core must remain inside its original context and contig: "+e.id);
		final int start=e.strand=='+'?e.originalStart-e.left+e.from-1:e.originalStop+e.left-e.to+1;
		final int stop=e.strand=='+'?e.originalStart-e.left+e.to-1:e.originalStop+e.left-e.from+1;
		require(start==e.start && stop==e.stop, "Sense-to-genome projection changed between the join and overlay: "+e.id);
	}
	private void write(String contexts, String out){
		RrnaResourceIO.fresh(out, ".overlay.tsv.gz", ".admitted.fa", ".followup.contexts.fa", ".followup.loci.tsv", ".summary.tsv");
		final ByteStreamWriter ledger=RrnaResourceIO.writer(out+".overlay.tsv.gz"), admitted=RrnaResourceIO.writer(out+".admitted.fa");
		final ByteStreamWriter followup=RrnaResourceIO.writer(out+".followup.contexts.fa"), loci=RrnaResourceIO.writer(out+".followup.loci.tsv");
		final ByteBuilder b=new ByteBuilder();final int[] counts=new int[5];
		try{
			ledger.println(HEADER);
			loci.println(NcrnaCmLocusAudit.HEADER);
			RrnaResourceIO.stream(contexts, r->{
				final Entry e=byId.get(RrnaResourceIO.first(r.id));require(e!=null && !e.seen, "Foreign or repeated source context: "+r.id);
				require(RrnaResourceIO.field(r.id, "source").equals(e.genome) && RrnaResourceIO.field(r.id, "contig").equals(e.seqid)
					&& RrnaResourceIO.field(r.id, "strand").equals(String.valueOf(e.strand))
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "start"))==e.originalStart
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "stop"))==e.originalStop
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "lflank"))==e.left
					&& Integer.parseInt(RrnaResourceIO.field(r.id, "rflank"))==e.right
					&& r.bases!=null && r.length()==e.length+e.left+e.right, "Context provenance or geometry changed: "+e.id);
				e.seen=true;counts[0]++;
				final boolean full=e.model>=0;final int length=full?e.to-e.from+1:0;
				if(full){
					require(length>=17, "A supported core shorter than the shared seed requires explicit policy: "+e.id);
					for(int i=e.from-1; i<e.to; i++){require(AminoAcid.baseToNumber[r.bases[i]&255]>=0, "Ambiguous supported core needs explicit disposition, not silent fitting: "+e.id);}
					counts[1]++;if(e.label.equals("CROSS_LABELED")){counts[2]++;}if(length<low || length>high){counts[4]++;}
					b.clear().append('>').append(e.id).append(" source_family=").append(family).append(" member_family=")
						.append(e.label.equals("CROSS_LABELED")?"B12clade":family).append(" source=").append(e.genome)
						.append(" contig=").append(e.seqid).append(" strand=").append(e.strand).append(" start=").append(e.start)
						.append(" stop=").append(e.stop).append(" coordinate_model=").append(MODELS[e.model]).append(" lflank=0 rflank=0\n");
					for(int i=e.from-1; i<e.to; i++){b.append((char)Character.toUpperCase((char)r.bases[i]));}admitted.print(b.nl());
				}else{
					counts[3]++;followup.print(b.clear().append('>').append(r.id).nl().append(r.bases).nl());loci.print(b.clear().append(e.source).nl());
				}
				b.clear().append(e.source).tab().append(family).tab().append(e.label).tab().append(full).tab().append(e.support)
					.tab().append(full?MODELS[e.model]:"NA").tab().append(e.stage).tab();
				if(full){b.append(e.start).tab().append(e.stop).tab().append(length).tab().append(length<low || length>high?"OUTSIDE_BAND_SUPPORTED":"WITHIN_BAND")
					.tab().append(e.label.equals("CROSS_LABELED")?"B12clade":family);}
				else{b.append("NA\tNA\tNA\tNA\tUNRESOLVED");}
				ledger.print(b.tab().append(!full?"UNRESOLVED":e.start==e.originalStart && e.stop==e.originalStop?"SOURCE_EXACT":"GEOMETRY-CORRECTED").nl());
			});
			require(counts[0]==entries.size() && counts[0]==counts[1]+counts[3], "Every source locus must be admitted or retained unresolved, exactly once");
			for(Entry e:entries){require(e.seen, "Original context missing from overlay: "+e.id);}
		}finally{RrnaResourceIO.finish(ledger, admitted, followup, loci);}
		final ByteStreamWriter summary=RrnaResourceIO.writer(out+".summary.tsv");
		try{
			summary.println("source_family\tpopulation\tadmitted\town_subtype\tcross_labeled\tunresolved\tadmitted_outside_length_band\tminlen\tmaxlen\tstatus");
			summary.print(b.clear().append(family).tab().append(counts[0]).tab().append(counts[1]).tab().append(counts[1]-counts[2])
				.tab().append(counts[2]).tab().append(counts[3]).tab().append(counts[4]).tab().append(low).tab().append(high).append("\tREVIEW_REQUIRED_BEFORE_FITTING\n"));
		}finally{RrnaResourceIO.finish(summary);}
		System.err.println("NCRNA_CLAN_OVERLAY_PASS population="+counts[0]+" admitted="+counts[1]+" unresolved="+counts[3]);
	}
	private static byte[] prefix(byte[] row, int fields){
		assert(fields>0):"Original source prefix must retain every annotation field";
		int n=0;for(int i=0; i<row.length; i++){if(row[i]=='\t' && ++n==fields){return Arrays.copyOf(row, i);}}
		throw new IllegalArgumentException("Incomplete source prefix");
	}
	private static String header(){
		final StringBuilder b=new StringBuilder(NcrnaCmLocusAudit.HEADER);
		b.append("\tsource_family\tfull_support_count\tsole_full_support_model\town_full_support\tevidence_class");
		for(String m:MODELS){for(String f:FIELDS){b.append('\t').append(m).append('_').append(f);}b.append('\t').append(m).append("_genomic_start1\t").append(m).append("_genomic_stop1");}
		return b.append("\tinterpretation").toString();
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Entry {
		Entry(String id_, byte[] source_, LineParser1 p){
			id=id_;source=source_;genome=p.parseString(1);seqid=p.parseString(6);length=p.parseInt(9);left=p.parseInt(18);right=p.parseInt(19);contigLength=p.parseInt(15);
			final String s=p.parseString(7), interval=p.parseString(8);final int dash=interval.indexOf('-');
			require((s.equals("+") || s.equals("-")) && dash>0 && interval.indexOf(',')<0, "Overlay requires explicit strand and contiguous source intervals");
			strand=s.charAt(0);originalStart=Integer.parseInt(interval.substring(0, dash));originalStop=Integer.parseInt(interval.substring(dash+1));
			require(originalStart>0 && originalStop<=contigLength && originalStop-originalStart+1==length && left>=0 && right>=0, "Invalid original context geometry");
		}
		final String id, genome, seqid;final byte[] source;final char strand;
		final int length, left, right, contigLength, originalStart, originalStop;
		int model=-1, from, to, start, stop;String label, support, stage;boolean initialUnresolved, rescueSeen, seen;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields/constants       ----------------*/
	/*--------------------------------------------------------------*/

	private final String family, phase;private final int own, low, high;private int unresolved;
	public static final String HEADER=NcrnaCmLocusAudit.HEADER+"\tsource_family\ttruth_subtype\tadmitted\tsupport_models\tcoordinate_model\tevidence_stage\tcurrent_start1\tcurrent_stop1\tcore_length\tdiagnostic_length_band\tmember_family\tgeometry_status";
	private final ArrayList<Entry> entries=new ArrayList<Entry>();
	private final ObjectMap<String,Entry> byId=new ObjectMap<String,Entry>(String.class, Entry.class);
	private static final String[] MODELS={"RF00174", "RF01482", "RF01689"};
	private static final int[] CLEN={189, 149, 81};private static final double[] GA={31.5, 50, 36.5};
	private static final String[] FIELDS={"cm_state", "raw_hits", "overlapping_sense_hits", "cm_model_from", "cm_model_to", "cm_seq_from", "cm_seq_to", "cm_trunc", "cm_bits", "overlap_bp", "left_excess", "right_excess", "geometry"};
}
