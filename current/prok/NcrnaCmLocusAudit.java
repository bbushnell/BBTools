package prok;

import java.util.ArrayList;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Attributes independent CM hits to sense-oriented annotated-locus contexts.
 * Preserves all annotations; full model span is a conservative evidence tier,
 * not biological endpoint truth or permission to train on an input sequence.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmLocusAudit {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaCmLocusAudit.class, false);
		String loci=null, hits=null, out=null, model=null, family=null;
		int expected=-1, clen=-1;double ga=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value: "+arg);
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("loci")){loci=value;}
			else if(key.equals("hits")){hits=value;}
			else if(key.equals("outprefix")){out=value;}
			else if(key.equals("modelid")){model=value;}
			else if(key.equals("family")){family=value;}
			else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}
			else if(key.equals("clen")){clen=Parse.parseIntKMG(value);}
			else if(key.equals("gabits")){ga=Double.parseDouble(value);}
			else{throw new IllegalArgumentException("Unknown parameter: "+key);}
		}
		require(loci!=null && hits!=null && out!=null && expected>0 && Double.isFinite(ga),
			"Locus attribution needs explicit source ledger, CM table, denominator and model criteria");
		final RrnaCmTruthBuilder.Criteria criteria=new RrnaCmTruthBuilder.Criteria(model, family, clen, clen, ga);
		Shared.setThreads(1);
		new NcrnaCmLocusAudit().run(loci, hits, out, expected, criteria);
		Shared.closeStream(pp.outstream);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(final String loci, final String hits, final String out, final int expected,
		final RrnaCmTruthBuilder.Criteria criteria){
		assert(criteria.minSpan==criteria.clen):"This diagnostic reports full model span without a fitted coverage cutoff";
		load(loci, expected);
		final ArrayList<RrnaCmTruthBuilder.Hit> raw=RrnaCmTruthBuilder.readHits(
			new RrnaCmTruthBuilder.Genome("corpus", "contexts", hits, null), criteria);
		for(RrnaCmTruthBuilder.Hit hit:raw){
			final Locus locus=byId.get(hit.seqid);
			require(locus!=null, "CM target is absent from the annotated context ledger: "+hit.seqid);
			require(Math.max(hit.seqFrom, hit.seqTo)<=locus.contextLength,
				"CM coordinates exceed the bound context: "+hit.seqid+" length="+locus.contextLength);
			locus.hits.add(hit);
		}
		RrnaResourceIO.fresh(out, ".loci.tsv.gz", ".summary.tsv", ".model_ends.tsv", ".lengths.tsv");
		final ByteStreamWriter ledger=RrnaResourceIO.writer(out+".loci.tsv.gz"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteStreamWriter ends=RrnaResourceIO.writer(out+".model_ends.tsv");
		final ByteBuilder b=new ByteBuilder();final long[] counts=new long[STATES.length];
		final long[][] spans=new long[criteria.clen+1][criteria.clen+1];long events=0;
		final LengthProfile[][] profiles=new LengthProfile[STATES.length][5];
		final LengthProfile[] reference=new LengthProfile[4];
		try{
			ledger.println(HEADER+"\tcm_state\traw_hits\toverlapping_sense_hits\tcm_model_from\tcm_model_to\tcm_seq_from\tcm_seq_to\tcm_trunc\tcm_bits\toverlap_bp\tleft_excess\tright_excess\tgeometry");
			for(Locus locus:ordered){
				RrnaCmTruthBuilder.Hit best=null;int overlapping=0, bestOverlap=0;
				for(RrnaCmTruthBuilder.Hit hit:locus.hits){
					events++;
					if(hit.strand!='+'){continue;}
					final int overlap=overlap(hit, locus);
					if(overlap<=0){continue;}
					overlapping++;
					if(best==null || better(hit, overlap, best, bestOverlap)){best=hit;bestOverlap=overlap;}
				}
				final int state=state(locus, best, overlapping, criteria);
				locus.best=best;
				final int geom=best==null?4:geometry(best, locus);
				if(profiles[state][geom]==null){profiles[state][geom]=new LengthProfile("ALL", STATES[state], geom);}
				profiles[state][geom].add(locus);
				if(state==5 && locus.sourceFull){
					if(reference[geom]==null){reference[geom]=new LengthProfile("SOURCE_FULL_CM_FULL", STATES[state], geom);}
					reference[geom].add(locus);
				}
				counts[state]++;
				b.clear().append(locus.row).tab().append(STATES[state]).tab().append(locus.hits.size()).tab().append(overlapping).tab();
				if(best==null){b.append("NA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA");}
				else{
					spans[best.mdlFrom][best.mdlTo]++;
					b.append(best.mdlFrom).tab().append(best.mdlTo).tab().append(best.seqFrom).tab().append(best.seqTo)
						.tab().append(best.trunc).tab().append(best.score).tab().append(bestOverlap)
						.tab().append(best.seqFrom-locus.left-1).tab().append(locus.left+locus.length-best.seqTo);
				}
				ledger.print(b.tab().append(best==null?"NONE":GEOMETRY[geometry(best, locus)]).nl());
			}
			require(events==raw.size(), "Every raw CM hit must survive attribution, including neighbors and opposite-strand hits");
			summary.println("state\tloci\traw_hit_events");long total=0;
			for(int i=0;i<counts.length;i++){total+=counts[i];summary.print(b.clear().append(STATES[i]).tab().append(counts[i]).append("\tNA\n"));}
			require(total==expected, "CM classification must conserve the full annotation denominator");
			summary.print(b.clear().append("TOTAL\t").append(total).tab().append(events).nl());
			ends.println("model_from\tmodel_to\tleft_deficit\tright_deficit\tselected_diagnostic_loci_including_ambiguous");
			for(int from=1;from<=criteria.clen;from++){
				for(int to=from;to<=criteria.clen;to++){
					if(spans[from][to]>0){ends.print(b.clear().append(from).tab().append(to).tab().append(from-1).tab().append(criteria.clen-to).tab().append(spans[from][to]).nl());}
				}
			}
		}finally{RrnaResourceIO.finish(ledger, summary, ends);}
		final ByteStreamWriter lengths=RrnaResourceIO.writer(out+".lengths.tsv");
		try{
			lengths.println("tier\tstate\tgeometry\tmetric\tcount\tmin\tp10\tmedian\tp90\tmax\tmin_id\tmin_left_excess\tmin_right_excess\tmax_id\tmax_left_excess\tmax_right_excess");
			for(LengthProfile[] group:profiles){for(LengthProfile p:group){if(p!=null){p.write(lengths, b);}}}
			for(LengthProfile p:reference){if(p!=null){p.write(lengths, b);}}
		}finally{RrnaResourceIO.finish(lengths);}
		System.err.println("NCRNA_CM_LOCUS_AUDIT_PASS loci="+expected+" rawHits="+raw.size()+"; no source annotation edits or eligibility selection");
	}

	private void load(final String path, final int expected){
		assert(ordered.isEmpty()):"A locus audit instance consumes exactly one immutable source ledger";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==21, "Annotated corpus ledger must retain its21 source/context fields");
				final String id=p.parseString(0);final int length=p.parseInt(9), left=p.parseInt(18), right=p.parseInt(19);
				require(length>0 && left>=0 && right>=0, "Context span must contain a nonempty original locus: "+id);
				final String intervals=p.parseString(8);
				require(intervals.indexOf(',')<0, "Multipart loci need interval-aware CM attribution; do not use concatenated length as genomic span: "+id);
				final int dash=intervals.indexOf('-');require(dash>0, "Missing inclusive genomic interval: "+id);
				require(Integer.parseInt(intervals.substring(dash+1))-Integer.parseInt(intervals.substring(0, dash))+1==length,
					"Contiguous context geometry must match the source interval length: "+id);
				final Locus locus=new Locus(id, row, length, left, right, p.termEquals("ANNOTATION_FULL_COPY", 14));
				require(!byId.contains(id), "Duplicate locus ID would merge independent annotation evidence: "+id);
				byId.put(id, locus);ordered.add(locus);
			}
		}finally{require(!in.close(), "Annotated locus ledger read failed: "+path);}
		require(ordered.size()==expected, "Locus ledger denominator changed: "+ordered.size()+" expected="+expected);
	}

	private static int overlap(final RrnaCmTruthBuilder.Hit h, final Locus l){
		assert(h.strand=='+'):"Original loci occupy left+1..left+length only on the sense-oriented context";
		return Math.max(0, Math.min(h.seqTo, l.left+l.length)-Math.max(h.seqFrom, l.left+1)+1);
	}

	private static boolean better(final RrnaCmTruthBuilder.Hit h, final int overlap,
		final RrnaCmTruthBuilder.Hit old, final int oldOverlap){
		assert(overlap>0 && oldOverlap>0):"Only overlapping sense hits participate in deterministic diagnostic ranking";
		if(overlap!=oldOverlap){return overlap>oldOverlap;}
		final int span=h.mdlTo-h.mdlFrom, oldSpan=old.mdlTo-old.mdlFrom;
		if(span!=oldSpan){return span>oldSpan;}
		final int score=Double.compare(Double.parseDouble(h.score), Double.parseDouble(old.score));
		if(score!=0){return score>0;}
		return h.line<old.line;
	}

	private static int state(final Locus l, final RrnaCmTruthBuilder.Hit h, final int count,
		final RrnaCmTruthBuilder.Criteria c){
		assert(count>=0 && (h==null)==(count==0)):"Presence of a diagnostic winner must agree with overlapping-hit count";
		if(h==null){return l.hits.isEmpty()?0:1;}
		if(count>1){return 2;}
		if(geometry(h, l)==3){return 6;}
		if(!h.included || Double.parseDouble(h.score)<c.ga){return 3;}
		if(!h.trunc.equals("no") || h.mdlFrom!=1 || h.mdlTo!=c.clen){return 4;}
		return 5;
	}

	private static int geometry(final RrnaCmTruthBuilder.Hit h, final Locus l){
		assert(h.strand=='+'):"Containment compares inclusive coordinates on the sense-oriented context";
		final boolean hitInside=h.seqFrom>=l.left+1 && h.seqTo<=l.left+l.length;
		final boolean annotationInside=h.seqFrom<=l.left+1 && h.seqTo>=l.left+l.length;
		return hitInside?(annotationInside?0:1):(annotationInside?2:3);
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	private static final class Locus {
		Locus(String id_, byte[] row_, int length_, int left_, int right_, boolean sourceFull_){
			id=id_;row=row_;length=length_;left=left_;contextLength=length_+left_+right_;sourceFull=sourceFull_;
		}
		final String id;final byte[] row;final int length, left, contextLength;final boolean sourceFull;
		RrnaCmTruthBuilder.Hit best;
		final ArrayList<RrnaCmTruthBuilder.Hit> hits=new ArrayList<RrnaCmTruthBuilder.Hit>();
	}

	/** Quantiles use sorted rank floor((n-1)*percent/100); extrema ties use locus ID. */
	private static final class LengthProfile {
		LengthProfile(String tier_, String state_, int geometry_){tier=tier_;state=state_;geom=geometry_;}
		void add(final Locus l){
			assert(l!=null && l.length>0):"Length summaries retain actual annotated loci, never synthetic lengths";
			source.add(l.length);srcMin=extreme(srcMin, l, false, true);srcMax=extreme(srcMax, l, false, false);
			if(l.best!=null){cm.add(l.best.seqTo-l.best.seqFrom+1);cmMin=extreme(cmMin, l, true, true);cmMax=extreme(cmMax, l, true, false);}
		}
		void write(final ByteStreamWriter out, final ByteBuilder b){
			assert(source.size>0):"Only observed state/geometry groups get length profiles";
			source.sort();cm.sort();writeMetric(out, b, "SOURCE", source, srcMin, srcMax);writeMetric(out, b, "CM_HIT", cm, cmMin, cmMax);
		}
		void writeMetric(final ByteStreamWriter out, final ByteBuilder b, final String metric,
			final IntList values, final Locus min, final Locus max){
			assert((values.size==0)==(min==null && max==null)):"Empty CM evidence must not invent extrema";
			b.clear().append(tier).tab().append(state).tab().append(geom==4?"NONE":GEOMETRY[geom]).tab().append(metric).tab().append(values.size);
			if(values.size==0){for(int i=0;i<5;i++){b.append("\tNA");}}
			else{for(int p:new int[]{0, 10, 50, 90, 100}){b.tab().append(values.array[(int)((values.size-1L)*p/100)]);}}
			appendExtreme(b, min);appendExtreme(b, max);out.print(b.nl());
		}
		static Locus extreme(final Locus old, final Locus next, final boolean hit, final boolean minimum){
			assert(next!=null && (!hit || next.best!=null)):"CM extrema require an observed selected hit";
			if(old==null){return next;}
			final int a=hit?next.best.seqTo-next.best.seqFrom+1:next.length, b=hit?old.best.seqTo-old.best.seqFrom+1:old.length;
			return a==b?(next.id.compareTo(old.id)<0?next:old):((minimum?a<b:a>b)?next:old);
		}
		static void appendExtreme(final ByteBuilder b, final Locus l){
			assert(b!=null):"Extrema and excess lengths share the same output row";
			if(l==null){b.append("\tNA\tNA\tNA");return;}
			b.tab().append(l.id);
			if(l.best==null){b.append("\tNA\tNA");}
			else{b.tab().append(l.best.seqFrom-l.left-1).tab().append(l.left+l.length-l.best.seqTo);}
		}
		final String tier, state;final int geom;
		final IntList source=new IntList(), cm=new IntList();
		Locus srcMin, srcMax, cmMin, cmMax;
	}

	private final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ArrayList<Locus> ordered=new ArrayList<Locus>();
	private static final String[] STATES={"NO_REPORTED_HIT", "NO_OVERLAPPING_SENSE_HIT", "MULTIPLE_OVERLAPPING_HITS",
		"BELOW_GA_OR_NOT_INCLUDED", "CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED", "CM_FULL_SPAN_SUPPORTED", "PARTIAL_ANNOTATION_OVERLAP_UNRESOLVED"};
	private static final String[] GEOMETRY={"EXACT", "HIT_WITHIN_ANNOTATION", "ANNOTATION_WITHIN_HIT", "PARTIAL_OVERLAP"};
	public static final String HEADER="locus_id\tgenome\tassembly\ttaxid\tdomain\tphylum\tseqid\tstrand\tintervals1\tlength\tsource_line\tsource_ids\tpartial\tpseudo\tstatus\tcontig_length\tedge\tundefined_bases\tlflank\trflank\tsource_lines";
}
