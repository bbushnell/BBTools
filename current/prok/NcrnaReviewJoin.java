package prok;

import java.util.ArrayList;
import java.util.Locale;
import java.util.function.Consumer;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins frozen prediction memberships with independent CM and source-feature evidence.
 * Evidence classes are not biological verdicts or permission to alter a truth set.
 * @author Brian, Raiden
 */
public final class NcrnaReviewJoin {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaReviewJoin.class, false);
		String members=null, context=null, overlaps=null, hits=null, cm=null, out=null;
		int expected=-1, events=-1, runs=-1, clen=-1;double ga=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("members")){members=value;}else if(key.equals("context")){context=value;}
			else if(key.equals("overlaps")){overlaps=value;}else if(key.equals("hits")){hits=value;}
			else if(key.equals("cm")){cm=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else if(key.equals("events")){events=Integer.parseInt(value);}
			else if(key.equals("runs")){runs=Integer.parseInt(value);}else if(key.equals("clen")){clen=Integer.parseInt(value);}
			else if(key.equals("gabits")){ga=Double.parseDouble(value);}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(members!=null && context!=null && overlaps!=null && hits!=null && cm!=null && out!=null
			&& expected>0 && events>=expected && runs>0 && clen>0 && Double.isFinite(ga), "Bind all five frozen tables and their declared dimensions/CM criteria");
		Shared.setThreads(1);
		try{
			final NcrnaReviewJoin x=new NcrnaReviewJoin(clen, ga);
			x.context(context, expected);x.features(overlaps);x.hits(hits);x.cm(cm);x.members(members, events, runs);x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaReviewJoin(int clen_, double ga_){clen=clen_;ga=ga_;}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void context(String path, int expected){
		read(path, CONTEXT_HEADER, p->{
			final Locus l=new Locus(p);
			require(!ids.contains(l.id) && l.modelLength==clen, "Unique review IDs and identical CM length are required: "+l.id);
			ids.put(l.id, l);loci.add(l);
		});
		require(loci.size()==expected, "Context denominator changed: "+loci.size()+" expected "+expected);
	}
	private void features(String path){
		read(path, OVERLAP_HEADER, p->{
			final Locus l=locus(p);final String kind=p.parseString(8), relation=p.parseString(11);
			require(p.termEquals(l.genome, 1) && p.termEquals(l.assembly, 2) && p.termEquals(l.fasta, 3) && p.termEquals(l.fastaHash, 4)
				&& p.termEquals(l.gff, 5) && p.termEquals(l.gffHash, 6) && p.termEquals(l.seqid, 13), "Feature source differs from frozen context: "+l.id);
			final int line=p.parseInt(7), overlap=p.parseInt(10), start=p.parseInt(16), stop=p.parseInt(17);
			require(line>0 && overlap>0 && overlap<=l.length && start>0 && stop>=start && l.featureLines.add(line), "Invalid/repeated source overlap row: "+l.id);
			final String strand=p.parseString(19);
			require(strand.equals("+") || strand.equals("-") || strand.equals(".") || strand.equals("?"), "Invalid GFF strand");
			final String want=strand.equals(".") || strand.equals("?")?"UNKNOWN":strand.equals(l.strand)?"SAME":"OPPOSITE";
			require(relation.equals(want), "Feature strand relation changed: "+l.id);
			if(kind.equals("REGION_METADATA")){l.seenRegions++;return;}
			l.seenFeatures++;
			if(kind.equals("PARENT_ANNOTATION")){l.seenParents++;return;}
			if(kind.equals("TARGET_RNA_ANNOTATION")){l.seenTarget++;}else if(kind.equals("OTHER_NCRNA")){l.seenOtherRna++;}
			else if(kind.equals("CDS")){l.seenCds++;}else{require(kind.equals("OTHER_FEATURE"), "Unknown source feature class: "+kind);l.seenOther++;}
			if(!relation.equals("SAME")){return;}
			final String attrs=p.parseString(21);final boolean target=rfam(attrs, "RF00174");
			if(target){l.targetSame++;}if(rfam(attrs, "RF01482")){l.variant1482++;}if(rfam(attrs, "RF01689")){l.variant1689++;}
			require(p.termEquals(target?"true":"false", 9), "RF00174 tag must match the exact Dbxref token");
		});
		for(Locus l:loci){
			require(l.seenFeatures==l.features && l.seenRegions==l.regions && l.seenParents==l.parents && l.seenTarget==l.target
				&& l.seenOtherRna==l.otherRna && l.seenCds==l.cds && l.seenOther==l.other, "Raw feature rows do not reproduce the frozen summary: "+l.id);
			require((l.targetSame>0?1:0)+(l.variant1482>0?1:0)+(l.variant1689>0?1:0)<=1, "Ambiguous same-strand B12 subtype annotations: "+l.id);
		}
	}
	private void hits(String path){
		read(path, HIT_HEADER, p->{
			final Locus l=locus(p);final int line=p.parseInt(1), start=p.parseInt(3), stop=p.parseInt(4), from=p.parseInt(6), to=p.parseInt(7);
			require(p.termEquals(l.seqid, 2) && line>0 && l.hitLines.add(line) && start>0 && stop>=start && from>=1 && to>=from && to<=clen,
				"Invalid/duplicate CM row or mismatched contig: "+l.id);
			final String strand=p.parseString(5);require(strand.equals("+") || strand.equals("-"), "CM strand must be explicit");
			final int overlap=strand.equals(l.strand)?Math.max(0, Math.min(stop, l.stop)-Math.max(start, l.start)+1):0;
			require(overlap==p.parseInt(12), "CM overlap must use exact call coordinates and strand: "+l.id);
			final double bits=p.parseDouble(9);require(Double.isFinite(bits), "CM score must be finite");
			final boolean included=bool(p, 11), full=included && bits>=ga && from==1 && to==clen && p.termEquals("no", 8) && overlap>0;
			final boolean strict=full && 2L*overlap>=l.length && 2L*overlap>=stop-start+1;
			require(bool(p, 13)==strict, "Stored CM qualification differs from declared GA/span/geometry: "+l.id);
			l.rawHits++;if(overlap>0){l.overlapHits++;}if(full){l.fullHits++;}if(strict){l.strictHits++;}
			if(overlap==0){return;}
			// Rank for display only; every raw hit remains in the bound input table.
			final int tier=strict?2:full?1:0;
			if(l.best==null || tier>l.bestTier || (tier==l.bestTier && bits>l.bestBits)){
				l.bestTier=tier;l.bestBits=bits;l.best=new ByteBuilder().append(start).tab().append(stop).tab().append(strand)
					.tab().append(from).tab().append(to).tab().append(p.parseString(8)).tab().append(bits, 3).tab().append(bits-ga, 3)
					.tab().append(overlap).tab().append(full).tab().append(strict).toBytes();
			}
		});
	}
	private void cm(String path){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		read(path, CM_HEADER, p->{
			final Locus l=locus(p);require(seen.add(l.id), "Duplicate CM summary: "+l.id);
			require(l.rawHits==p.parseInt(1) && l.overlapHits==p.parseInt(2) && l.strictHits==p.parseInt(3), "Raw hits must reproduce the frozen CM counts: "+l.id);
			final String status=l.strictHits==0?"NO_FULL_GA_OVERLAPPING_HIT":l.strictHits==1?"ONE_FULL_GA_OVERLAPPING_HIT":"MULTIPLE_FULL_GA_OVERLAPPING_HITS";
			require(p.termEquals(status, 4), "CM summary status disagrees with hit count: "+l.id);
		});
		require(seen.size()==loci.size(), "CM summary must contain every context locus");
	}
	private void members(String path, int expected, int runCount){
		final ObjectSet<String> seen=new ObjectSet<String>(String.class), names=new ObjectSet<String>(String.class);
		read(path, MEMBER_HEADER, p->{
			final String run=p.parseString(0), call=p.parseString(1), id=p.parseString(2);
			final Locus l=ids.get(id);require(l!=null && !run.isEmpty() && !call.isEmpty() && seen.add(run+'\t'+call), "Foreign or repeated membership: "+run+" "+call);
			final String source=p.parseString(3), full=p.parseString(4);
			require(outcome(source) && outcome(full) && (source.equals("APPARENT_FP") || full.equals("APPARENT_FP")), "Review membership must retain its original FP disposition");
			members.add(new Member(run, call, source, full, l));l.members++;names.add(run);
		});
		require(members.size()==expected && names.size()==runCount, "Review membership/run denominator changed");
		for(Locus l:loci){require(l.members>0, "Unreferenced context would silently broaden the union: "+l.id);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(String prefix){
		RrnaResourceIO.fresh(prefix, ".calls.tsv", ".summary.tsv");
		final ByteStreamWriter calls=RrnaResourceIO.writer(prefix+".calls.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ObjectMap<String,Counts> counts=new ObjectMap<String,Counts>(String.class, Counts.class);final ArrayList<Counts> ordered=new ArrayList<Counts>();
		final ByteBuilder b=new ByteBuilder();
		try{
			calls.println("run\tcall_id\tsource_outcome\tfull_outcome\t"+CONTEXT_HEADER+"\tsame_rf00174_rows\tsame_rf01482_rows\tsame_rf01689_rows\traw_cm_hits\tsense_cm_hits\tfull_ga_cm_hits\tstrict_cm_hits\tbest_cm_start1\tbest_cm_stop1\tbest_cm_strand\tbest_model_from\tbest_model_to\tbest_trunc\tbest_bits\tbest_ga_margin\tbest_overlap_bp\tbest_full_ga\tbest_strict\tevidence_class\tinterpretation");
			for(Member m:members){
				final Locus l=m.locus;final String category=l.category();
				calls.print(b.clear().append(m.run).tab().append(m.call).tab().append(m.source).tab().append(m.full).tab().append(l.raw)
					.tab().append(l.targetSame).tab().append(l.variant1482).tab().append(l.variant1689).tab().append(l.rawHits).tab().append(l.overlapHits)
					.tab().append(l.fullHits).tab().append(l.strictHits).tab().append(l.best==null?NO_HIT:l.best).tab().append(category)
					.append("\tEVIDENCE_ONLY_NO_TRUTH_EDIT\n"));
				final String key=m.run+'\t'+category;Counts c=counts.get(key);
				if(c==null){c=new Counts(m.run, category);counts.put(key, c);ordered.add(c);}c.events++;
				if(m.source.equals("APPARENT_FP")){c.sourceFp++;}if(m.full.equals("APPARENT_FP")){c.fullFp++;}c.loci.add(l.id);
			}
			summary.println("run\tevidence_class\tmemberships\tdistinct_exact_loci\tsource_FP\tfull_FP\tinterpretation");
			int total=0;for(Counts c:ordered){
				total+=c.events;summary.print(b.clear().append(c.run).tab().append(c.category).tab().append(c.events).tab().append(c.loci.size())
					.tab().append(c.sourceFp).tab().append(c.fullFp).append("\tNOT_BIOLOGICAL_TP_OR_FP_VERDICTS\n"));
			}
			require(total==members.size(), "Evidence partition must conserve every membership, without deduplicating biological loci by guesswork");
		}finally{RrnaResourceIO.finish(calls, summary);}
		System.err.println("NCRNA_REVIEW_JOIN_PASS loci="+loci.size()+" memberships="+members.size());
	}

	/*--------------------------------------------------------------*/
	/*----------------         Static methods       ----------------*/
	/*--------------------------------------------------------------*/

	private Locus locus(LineParser1 p){final Locus l=ids.get(p.parseString(0));require(l!=null, "Unknown review ID in evidence");return l;}
	private static void read(String path, String header, Consumer<LineParser1> action){
		assert(path!=null && header!=null):"Every join input has a declared source and exact schema";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');p.set(header.getBytes(java.nio.charset.StandardCharsets.UTF_8));final int terms=p.terms();
		try{RrnaResourceIO.header(in, header, path);for(byte[] row; (row=in.nextLine())!=null; ){p.set(row);require(p.terms()==terms, "Wrong field count: "+path);action.accept(p);}}
		finally{require(!in.close(), "Failed reading join input: "+path);}
	}
	private static boolean rfam(String attrs, String model){
		assert(attrs!=null && model!=null):"Subtype evidence uses exact Rfam Dbxref tokens";
		final String value=RrnaResourceIO.optionalField(attrs, "Dbxref");if(value==null){return false;}
		for(int from=0; from<value.length(); ){
			int to=value.indexOf(',', from);if(to<0){to=value.length();}
			if(to-from==5+model.length() && value.regionMatches(true, from, "Rfam:", 0, 5) && value.regionMatches(from+5, model, 0, model.length())){return true;}from=to+1;
		}return false;
	}
	private static boolean outcome(String s){return s.equals("TP") || s.equals("APPARENT_FP") || s.equals("UNSCORED_UNRESOLVED_OVERLAP");}
	private static boolean bool(LineParser1 p, int col){require(p.termEquals("true", col) || p.termEquals("false", col), "Malformed boolean in frozen evidence");return p.termEquals("true", col);}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Locus {
		Locus(LineParser1 p){
			id=p.parseString(0);genome=p.parseString(1);assembly=p.parseString(2);seqid=p.parseString(3);start=p.parseInt(4);stop=p.parseInt(5);strand=p.parseString(6);
			length=p.parseInt(7);modelLength=p.parseInt(8);features=p.parseInt(11);regions=p.parseInt(12);parents=p.parseInt(13);target=p.parseInt(14);
			otherRna=p.parseInt(15);cds=p.parseInt(16);other=p.parseInt(17);context=p.parseString(21);gff=p.parseString(22);gffHash=p.parseString(23);fasta=p.parseString(24);fastaHash=p.parseString(25);
			require(start>0 && stop>=start && length==stop-start+1 && modelLength>0 && (strand.equals("+") || strand.equals("-"))
				&& features>=0 && regions>=0 && parents>=0 && target>=0 && otherRna>=0 && cds>=0 && other>=0
				&& features==parents+target+otherRna+cds+other && features==p.parseInt(18)+p.parseInt(19)+p.parseInt(20)
				&& p.parseInt(9)==length-modelLength && Math.abs(p.parseDouble(10)-(double)length/modelLength)<0.000001,
				"Context geometry/type/strand accounting must be internally consistent: "+id);
			require(gffHash.matches("[0-9a-f]{20}") && fastaHash.matches("[0-9a-f]{20}"), "Context provenance must be sha80");
			require(context.equals(features==0?"INTERGENIC":features==parents?"PARENT_ONLY_OR_AMBIGUOUS":"FEATURE_OVERLAP"), "Context label disagrees with feature counts");
			final ByteBuilder b=new ByteBuilder();for(int i=0; i<p.terms(); i++){if(i>0){b.tab();}b.append(p.parseByteArray(i));}raw=b.toBytes();
		}
		String category(){
			assert(members>0):"Only original review memberships may enter the evidence report";
			if(variant1482+variant1689>0){return "OTHER_B12_ANNOTATION_SCOPE_UNRESOLVED";}
			if(targetSame>0){return "TARGET_ANNOTATION_GEOMETRY_REVIEW";}
			if(strictHits>0){return features==0?"FULL_CM_INTERGENIC_SUPPORT":"FULL_CM_WITH_FEATURE_OVERLAP";}
			if(fullHits>0){return "FULL_CM_GEOMETRY_REVIEW";}
			if(overlapHits>0){return "PARTIAL_OR_BELOW_GA_CM_AMBIGUOUS";}
			return features==0?"NO_REPORTED_CM_INTERGENIC_AMBIGUOUS":"NO_REPORTED_CM_FEATURE_OVERLAP_AMBIGUOUS";
		}
		final String id, genome, assembly, seqid, strand, context, gff, gffHash, fasta, fastaHash;
		final int start, stop, length, modelLength, features, regions, parents, target, otherRna, cds, other;
		final byte[] raw;byte[] best;double bestBits;int bestTier;
		int seenFeatures, seenRegions, seenParents, seenTarget, seenOtherRna, seenCds, seenOther, targetSame, variant1482, variant1689;
		int rawHits, overlapHits, fullHits, strictHits, members;
		final map.IntHashSet featureLines=new map.IntHashSet(), hitLines=new map.IntHashSet();
	}
	private static final class Member {
		Member(String r, String c, String s, String f, Locus l){run=r;call=c;source=s;full=f;locus=l;}
		final String run, call, source, full;final Locus locus;
	}
	private static final class Counts {
		Counts(String r, String c){run=r;category=c;}
		final String run, category;int events, sourceFp, fullFp;final ObjectSet<String> loci=new ObjectSet<String>(String.class);
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields/constants       ----------------*/
	/*--------------------------------------------------------------*/

	private final int clen;private final double ga;
	private final ObjectMap<String,Locus> ids=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<Member> members=new ArrayList<Member>();
	private static final byte[] NO_HIT="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA".getBytes(java.nio.charset.StandardCharsets.US_ASCII);
	public static final String MEMBER_HEADER="run\tcall_id\treview_id\tsource_outcome\tfull_outcome";
	public static final String CM_HEADER="review_id\traw_hits\tsense_overlapping_hits\tfull_ga_reciprocal_half_hits\tstatus";
	public static final String HIT_HEADER="review_id\tcm_line\tseqid\tgenomic_start1\tgenomic_stop1\tgenomic_strand\tmodel_from\tmodel_to\ttrunc\tbits\tevalue\tincluded\tsense_overlap_bp\tfull_ga_reciprocal_half";
	public static final String CONTEXT_HEADER="review_id\tgenome\tassembly\tseqid\tstart1\tstop1\tstrand\tcall_length\tmodel_length\tlength_minus_model\tlength_ratio\tfeature_rows\tregion_rows\tparent_rows\ttarget_annotation_rows\tother_ncrna_rows\tcds_rows\tother_rows\tsame_strand_rows\topposite_strand_rows\tunknown_strand_rows\tcontext\tgff\tgff_sha80\tfasta\tfasta_sha80";
	public static final String OVERLAP_HEADER="review_id\tgenome\tassembly\tfasta\tfasta_sha80\tgff\tgff_sha80\tgff_line\tfeature_class\ttarget_rfam_tag\toverlap_bp\tstrand_relation\torigin_spanning\tseqid\tsource\ttype\tstart1\tstop1\tscore\tstrand\tphase\toriginal_attributes";
}
