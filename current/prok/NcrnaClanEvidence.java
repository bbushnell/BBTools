package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;

/** Joins model-specific CM audits without assigning an ambiguous biological subtype.
 * Legacy audits retain their single-hit containment tier. Raw adjudication also
 * supports explicit any-hit/full-model evidence without requiring containment.
 * Proposed bounds reproduce the observed own-model/exact/source-full extrema rule.
 * @author Brian, Raiden
 */
public final class NcrnaClanEvidence {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaClanEvidence.class, false);
		String manifest=null, family=null, out=null;int expected=-1, modelCount=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("models")){manifest=value;}else if(key.equals("sourcefamily")){family=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else if(key.equals("modelcount")){modelCount=Integer.parseInt(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(manifest!=null && family!=null && family.matches("RF[0-9]{5}") && out!=null && expected>0 && modelCount>0,
			"Bind the model/audit manifest, source Rfam family and complete locus/model denominators");
		Shared.setThreads(1);
		try{
			final NcrnaClanEvidence x=new NcrnaClanEvidence(family);
			x.models(manifest, modelCount);for(int m=0; m<x.models.size(); m++){x.audit(m, expected);}x.write(out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaClanEvidence(String family_){family=family_;}
	private void models(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> names=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, "model_id\tclen\tga_bits\taudit", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==4, "Each model needs its own criteria and complete audit");
				final String id=p.parseString(0), audit=p.parseString(3);final int clen=p.parseInt(1);final double ga=p.parseDouble(2);
				require(id.matches("RF[0-9]{5}") && names.add(id) && clen>0 && Double.isFinite(ga) && audit.startsWith("/"), "Invalid or duplicated model binding");
				if(id.equals(family)){own=models.size();}models.add(new Model(id, audit, clen, ga));
			}
		}finally{require(!in.close(), "Model manifest read failed");}
		require(models.size()==expected && own>=0, "Every declared model, including the source family's own CM, must be present");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void audit(int index, int expected){
		assert(index>=0 && index<models.size()):"Every audit occupies exactly one declared model column";
		final Model model=models.get(index);final ByteFile in=RrnaResourceIO.open(model.path);final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			RrnaResourceIO.header(in, NcrnaCoreRecovery.AUDIT_HEADER, model.path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==34, "A clan join requires intact21source+13audit fields");
				final String id=p.parseString(0);require(id.startsWith(family+"_"), "Audit source family changed: "+id);
				final int boundary=sourceEnd(row);final byte[] source=Arrays.copyOf(row, boundary);
				Entry e=byId.get(id);
				if(index==0){require(e==null, "Duplicate source ID: "+id);e=new Entry(id, source, p, models.size());byId.put(id, e);entries.add(e);}
				require(e!=null && e.audit[index]==null && Arrays.equals(e.source, source), "Missing/repeated ID or changed source metadata across model audits: "+id);
				final String state=p.parseString(21);require(validState(state), "Unknown CM evidence state: "+state);
				e.audit[index]=Arrays.copyOfRange(row, boundary+1, row.length);count++;
				final boolean any=p.termEquals("CM_FULL_SPAN_ANY_HIT_SUPPORTED", 21);
				final boolean full=any || p.termEquals("CM_FULL_SPAN_SUPPORTED", 21);e.full[index]=full;
				if(full){
					require((any?p.parseInt(23)>=1:p.parseInt(23)==1) && p.parseInt(24)==1 && p.parseInt(25)==model.clen && p.termEquals("no", 28) && p.parseDouble(29)>=model.ga && p.parseInt(30)>0,
						"Full-support label must retain an overlapping untruncated GA-passing complete-model hit: "+id);
					require(p.termEquals("EXACT", 33) || p.termEquals("HIT_WITHIN_ANNOTATION", 33) || p.termEquals("ANNOTATION_WITHIN_HIT", 33)
						|| (any && p.termEquals("PARTIAL_OVERLAP", 33)), "Only explicit raw adjudication retires the legacy containment tier");
				}
				if(!p.termEquals("NA", 26)){
					final int a=p.parseInt(26), z=p.parseInt(27);
					require(a>0 && z>=a && z<=e.length+e.left+e.right, "Selected CM interval is outside its original sense context: "+id);
					e.start[index]=e.strand=='+'?e.originalStart-e.left+a-1:e.originalStop+e.left-z+1;
					e.stop[index]=e.strand=='+'?e.originalStart-e.left+z-1:e.originalStop+e.left-a+1;
					require(e.start[index]>0 && e.stop[index]<=e.contigLength && e.stop[index]-e.start[index]==z-a, "CM projection must preserve length inside its source contig: "+id);
				}else{require(!full, "A full CM hit cannot lack sequence coordinates");}
				if(index==own && full && e.sourceEligible && p.termEquals("EXACT", 33)){
					lengths.add(e.length);
					if(min==null || e.length<min.length || (e.length==min.length && id.compareTo(min.id)<0)){min=e;}
					if(max==null || e.length>max.length || (e.length==max.length && id.compareTo(max.id)<0)){max=e;}
				}
			}
		}finally{require(!in.close(), "CM audit read failed: "+model.path);}
		require(count==expected && entries.size()==expected, "Every model must cover the same complete source population");
		for(Entry e:entries){require(e.audit[index]!=null, "Model audit silently omitted a locus: "+e.id);}
	}
	private void write(String prefix){
		RrnaResourceIO.fresh(prefix, ".loci.tsv.gz", ".summary.tsv", ".bounds.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(prefix+".loci.tsv.gz"), summary=RrnaResourceIO.writer(prefix+".summary.tsv"), bounds=RrnaResourceIO.writer(prefix+".bounds.tsv");
		final ByteBuilder b=new ByteBuilder();final int[] counts=new int[CATEGORIES.length];
		try{
			b.append(NcrnaCmLocusAudit.HEADER).append("\tsource_family\tfull_support_count\tsole_full_support_model\town_full_support\tevidence_class");
			for(Model model:models){for(String field:AUDIT_FIELDS){b.tab().append(model.id).append('_').append(field);}b.tab().append(model.id).append("_genomic_start1").tab().append(model.id).append("_genomic_stop1");}
			rows.print(b.append("\tinterpretation\n"));
			for(Entry e:entries){
				int full=0, sole=-1;for(int m=0; m<models.size(); m++){require(e.audit[m]!=null, "Output cannot conceal an incomplete model join");if(e.full[m]){full++;sole=m;}}
				final int category=full==0?0:e.full[own]?(full==1?1:2):(full==1?3:4);counts[category]++;
				b.clear().append(e.source).tab().append(family).tab().append(full).tab().append(full==1?models.get(sole).id:"NA").tab().append(e.full[own]).tab().append(CATEGORIES[category]);
				for(int m=0; m<models.size(); m++){
					b.tab().append(e.audit[m]).tab();if(e.start[m]>0){b.append(e.start[m]).tab().append(e.stop[m]);}else{b.append("NA\tNA");}
				}
				rows.print(b.append("\tEVIDENCE_ONLY_SUBTYPE_AND_ADMISSION_REVIEW_REQUIRED\n"));
			}
			summary.println("source_family\tevidence_class\tloci\tinterpretation");int total=0;
			for(int i=0; i<counts.length; i++){total+=counts[i];summary.print(b.clear().append(family).tab().append(CATEGORIES[i]).tab().append(counts[i]).append("\tNOT_AUTOMATIC_ADMISSION_OR_REASSIGNMENT\n"));}
			require(total==entries.size(), "Support-pattern classes must conserve every original annotation");
			bounds.println("source_family\tclen\texact_full_source_loci\tproposed_min\tp10\tmedian\tp90\tproposed_max\tmin_id\tmax_id\tstatus");
			b.clear().append(family).tab().append(models.get(own).clen).tab().append(lengths.size);
			if(lengths.size==0){b.append("\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNO_EXACT_FULL_SOURCE_EVIDENCE\n");}
			else{
				lengths.sort();final int n=lengths.size;
				b.tab().append(min.length).tab().append(lengths.get((n-1)/10)).tab().append(lengths.get((n-1)/2))
					.tab().append(lengths.get((int)((n-1)*9L/10))).tab().append(max.length).tab().append(min.id).tab().append(max.id)
					.append("\tPROPOSED_EXTREMA_REVIEW_NOT_A_CALLER_CUTOFF\n");
			}
			bounds.print(b);
		}finally{RrnaResourceIO.finish(rows, summary, bounds);}
		System.err.println("NCRNA_CLAN_EVIDENCE_PASS source="+family+" loci="+entries.size()+" models="+models.size()+" own_exact="+lengths.size);
	}
	private static int sourceEnd(byte[] row){
		int fields=0;for(int i=0; i<row.length; i++){if(row[i]=='\t' && ++fields==21){return i;}}
		throw new IllegalArgumentException("Missing21-field source prefix");
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static boolean validState(String s){
		for(String state:AUDIT_STATES){if(s.equals(state)){return true;}}return false;
	}
	private static boolean flag(LineParser1 p, int column){require(p.termEquals("true", column) || p.termEquals("false", column), "Source flags must be explicit booleans");return p.termEquals("true", column);}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Model {
		Model(String i, String p, int c, double g){id=i;path=p;clen=c;ga=g;}
		final String id, path;final int clen;final double ga;
	}
	private static final class Entry {
		Entry(String id_, byte[] source_, LineParser1 p, int models){
			id=id_;source=source_;length=p.parseInt(9);left=p.parseInt(18);right=p.parseInt(19);contigLength=p.parseInt(15);
			final String interval=p.parseString(8), s=p.parseString(7);final int dash=interval.indexOf('-');
			require(dash>0 && interval.indexOf(',')<0 && (s.equals("+") || s.equals("-")), "A CM audit must retain contiguous source coordinates and explicit strand");
			originalStart=Integer.parseInt(interval.substring(0, dash));originalStop=Integer.parseInt(interval.substring(dash+1));strand=s.charAt(0);
			require(length>0 && originalStart>0 && originalStop<=contigLength && originalStop-originalStart+1==length && left>=0 && right>=0, "Source geometry is inconsistent: "+id);
			final boolean partial=flag(p, 12), pseudo=flag(p, 13), edge=flag(p, 16), sourceFull=p.termEquals("ANNOTATION_FULL_COPY", 14);
			final int undefined=p.parseInt(17);
			require(undefined>=0 && undefined<=length && (!sourceFull || (!partial && !pseudo)), "Source-full flags and ambiguity count must agree: "+id);
			require(sourceFull || p.termEquals("SOURCE_PARTIAL", 14) || p.termEquals("PSEUDOGENE", 14) || p.termEquals("MULTIPART_UNRESOLVED", 14), "Unknown source-copy status");
			sourceEligible=sourceFull && !edge && undefined==0;
			audit=new byte[models][];full=new boolean[models];start=new int[models];stop=new int[models];
		}
		final String id;final byte[] source;final int length, left, right, contigLength, originalStart, originalStop;final char strand;final boolean sourceEligible;
		final byte[][] audit;final boolean[] full;final int[] start, stop;
	}

	/*--------------------------------------------------------------*/
	/*----------------       Fields/constants       ----------------*/
	/*--------------------------------------------------------------*/

	private final String family;private int own=-1;
	private final ArrayList<Model> models=new ArrayList<Model>();private final ArrayList<Entry> entries=new ArrayList<Entry>();
	private final ObjectMap<String,Entry> byId=new ObjectMap<String,Entry>(String.class, Entry.class);
	private final IntList lengths=new IntList();private Entry min, max;
	private static final String[] CATEGORIES={"NO_UNAMBIGUOUS_FULL_SUPPORT", "OWN_FULL_ONLY", "OWN_AND_OTHER_FULL", "OTHER_FULL_ONLY", "MULTIPLE_OTHER_FULL"};
	private static final String[] AUDIT_STATES={"NO_REPORTED_HIT", "NO_OVERLAPPING_SENSE_HIT", "MULTIPLE_OVERLAPPING_HITS", "BELOW_GA_OR_NOT_INCLUDED", "CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED", "CM_FULL_SPAN_SUPPORTED", "PARTIAL_ANNOTATION_OVERLAP_UNRESOLVED", "CM_FULL_SPAN_ANY_HIT_SUPPORTED", "RAW_NO_QUALIFYING_FULL_HIT"};
	private static final String[] AUDIT_FIELDS={"cm_state", "raw_hits", "overlapping_sense_hits", "cm_model_from", "cm_model_to", "cm_seq_from", "cm_seq_to", "cm_trunc", "cm_bits", "overlap_bp", "left_excess", "right_excess", "geometry"};
}
