package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import idaligner.ScrabbleAligner;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.FloatList;

/** Derives model cutoffs from saved assigned-member identities, never best matches.
 * Explicit singleton contributors are checked against the original sequence.
 * Saved identities have six decimals; this does not reconstruct extra precision.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaMemberIdentityCutoffs {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaMemberIdentityCutoffs.class, false);
		String input=null, refs=null, quality=null, membership=null, additions=null, out=null, family=null;
		int expected=-1, expectedModels=-1;float floor=.60f;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){input=v;}else if(k.equals("ref")){refs=v;}else if(k.equals("quality")){quality=v;}
			else if(k.equals("membership")){membership=v;}else if(k.equals("additions")){additions=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("family")){family=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("models")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("floor")){floor=Float.parseFloat(v);}
			else{throw new IllegalArgumentException("Unknown parameter: "+k);}
		}
		require(input!=null && refs!=null && quality!=null && membership!=null && out!=null
			&& family!=null && family.matches("RF[0-9]{5}") && expected>0 && expectedModels>0,
			"Bind complete original corpus, final library, assigned identity report, membership report, family and both counts");
		require(Float.isFinite(floor) && floor>=0 && floor<=.80f, "Identity floor must be finite in [0,.80], at or below the fixed cap");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		final NcrnaMemberIdentityCutoffs exporter=new NcrnaMemberIdentityCutoffs();exporter.floor=floor;
		try{exporter.run(input, refs, quality, membership, additions, out, family, expected, expectedModels);}
		finally{Shared.closeStream(pp.outstream);}
	}

	private void run(String input, String refs, String quality, String membership, String additions,
		String out, String family, int expected, int expectedModels){
		RrnaResourceIO.stream(refs, r->{
			final String id=RrnaResourceIO.first(r.id);
			require(id.startsWith("cobalamin_"+family+'_') && !ids.contains(id) && r.bases!=null && r.length()>0,
				"Unique nonempty consensus names must belong to declared family: "+id);
			ids.put(id, models.size());models.add(new Model(r));
		});
		require(models.size()==expectedModels, "Final library count changed");
		readQuality(quality, expected);
		if(additions!=null){readAdditions(additions, family);}
		final int[] seen={0};
		RrnaResourceIO.stream(input, r->{
			final String id=RrnaResourceIO.first(r.id);final Member m=members.get(id);
			require(m!=null && !m.seen && r.length()==m.length, "Quality rows must bind exactly once to original input: "+id);
			m.seen=true;seen[0]++;
			if(m.added>=0){
				final Model model=models.get(m.added);
				require(Arrays.equals(r.bases, model.ref.bases) && RrnaResourceIO.field(r.id, "member_family").equals(m.addedLabel),
					"Singleton identity requires exact original sequence and biological label: "+id);
				// Tracer.postprocess gives literal N half credit even for byte-identical sequences.
				final float identity=ScrabbleAligner.alignStatic(r.bases, model.ref.bases, null);
				require(Float.isFinite(identity) && identity>=0 && identity<=1, "Singleton alignment must return a finite identity: "+id);
				model.add(m.length, identity);
			}
		});
		require(seen[0]==expected, "Missing original sequences cannot support construction identities");
		validateMembership(membership, family);
		write(out, family);
	}

	private void readQuality(String path, int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, QUALITY_HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==8, "Saved assigned identity schema changed");
				final String id=p.parseString(0);final int length=p.parseInt(1);
				require(!members.contains(id) && length>0, "Each quality member needs one nonempty sequence: "+id);
				final Member m=new Member(length);members.put(id, m);
				if(p.termEquals("kept", 2)){
					final String name=p.parseString(3);require(ids.contains(name), "Assigned model missing: "+name);
					final float identity=p.parseFloat(4);require(Float.isFinite(identity) && identity>=0 && identity<=1,
						"Assigned identity must be a finite fraction, not a missing or best-model value: "+id);
					models.get(ids.get(name)).add(length, identity);
				}else{
					require((p.termEquals("small", 2) || p.termEquals("orphan", 2)) && p.termEquals("NA", 3) && p.termEquals("NA", 4),
						"Unassigned members must not acquire a best-match construction assignment: "+id);
				}
			}
		}finally{require(!in.close(), "Assigned identity report read failed");}
		require(members.size()==expected, "Saved report must conserve all admitted members");
	}

	private void readAdditions(String path, String family){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, ADDITIONS_HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==6 && p.termEquals(family, 2) && p.termEquals("EXACT_SINGLE_MEMBER_NO_REASSIGNMENT", 5),
					"Singleton contributions need explicit original construction provenance");
				final Member m=members.get(p.parseString(0));final String name=p.parseString(1);
				require(m!=null && m.added<0 && ids.contains(name) && m.length==p.parseInt(4), "Invalid or duplicate singleton member: "+name);
				final int index=ids.get(name);final Model model=models.get(index);
				require(model.values.size==0 && !model.singleton && model.ref.length()==m.length,
					"A singleton cannot overwrite an original construction model: "+name);
				m.added=index;m.addedLabel=p.parseString(3);model.singleton=true;
			}
		}finally{require(!in.close(), "Singleton construction record read failed");}
	}

	private void validateMembership(String path, String family){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int rows=0;
		try{
			RrnaResourceIO.header(in, NcrnaMemberLengthBounds.HEADER, path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==19 && p.termEquals(family, 1) && ids.contains(p.parseString(0)), "Membership must name the exact final model set");
				final Model m=models.get(ids.get(p.parseString(0)));
				require(!m.checked && m.values.size>0 && m.ref.length()==p.parseInt(2) && m.values.size==p.parseInt(3)
					&& m.min==p.parseInt(4) && m.max==p.parseInt(8), "Assigned identities must conserve construction count and raw extrema: "+m.ref.id);
				m.checked=true;rows++;
			}
		}finally{require(!in.close(), "Membership cross-check read failed");}
		require(rows==models.size(), "Every model must have nonempty validated construction memberships");
	}

	private void write(String prefix, String family){
		RrnaResourceIO.fresh(prefix, "A.tsv", "B.tsv", "C.tsv", "D.tsv", "models.tsv", "distributions.tsv");
		final ByteStreamWriter detail=RrnaResourceIO.writer(prefix+"models.tsv"), summary=RrnaResourceIO.writer(prefix+"distributions.tsv");
		final ByteBuilder b=new ByteBuilder();final FloatList[] cuts={new FloatList(), new FloatList(), new FloatList(), new FloatList()};
		try{
			detail.println("model\tmembers\tsingleton\traw_min\traw_p5\traw_median\tA\tB\tC\tD\tquantile_definition\tidentity_source");
			for(Model m:models){
				m.values.sort();final float min=quantile(m.values, 0), p5=quantile(m.values, 5), median=quantile(m.values, 50);
				cuts[0].add(clamp(min));cuts[1].add(clamp(p5));cuts[2].add(clamp(p5-.02f));cuts[3].add(clamp(median-.05f));
				b.clear().append(RrnaResourceIO.first(m.ref.id)).tab().append(m.values.size).tab().append(m.singleton).tab().append(min, 6).tab().append(p5, 6).tab().append(median, 6);
				for(FloatList values:cuts){b.tab().append(values.array[values.size-1], 6);}
				detail.print(b.append("\tlower_order_statistic_floor((n-1)*p)\tassigned_Scrabble_shorter_query_saved_6_decimals_plus_exact_singletons\n"));
			}
			summary.println("family\trule\tmodels\tmin\tp5\tmedian\tp95\tmax\tat_lower_clamp\tat_upper_clamp\tfloor\tcap");
			for(int rule=0;rule<cuts.length;rule++){
				final String name=String.valueOf((char)('A'+rule));final FloatList values=cuts[rule];
				final ByteStreamWriter table=RrnaResourceIO.writer(prefix+name+".tsv");int lo=0, hi=0;
				try{
					table.println("model\tmodel_length\tmin_length\tmax_length\tidentity");
					for(int i=0;i<models.size();i++){
						final Model m=models.get(i);final float value=values.array[i];if(value==floor){lo++;}if(value==.80f){hi++;}
						// loadLibrary preserves r.id, including construction metadata after the first token.
						table.print(b.clear().append(m.ref.id).tab().append(m.ref.length())
							.tab().append(NcrnaMitoLengthProfile.lowerBound(m.min)).tab().append(NcrnaMitoLengthProfile.upperBound(m.max)).tab().append(value, 6).nl());
					}
				}finally{RrnaResourceIO.finish(table);}
				values.sort();b.clear().append(family).tab().append(name).tab().append(models.size());
				for(int q:new int[]{0, 5, 50, 95, 100}){b.tab().append(quantile(values, q), 6);}
				summary.print(b.tab().append(lo).tab().append(hi).tab().append(floor, 6).append("\t0.800000\n"));
			}
		}finally{RrnaResourceIO.finish(detail, summary);}
		System.err.println("NCRNA_MEMBER_IDENTITY_CUTOFFS_PASS models="+models.size()+" members="+members.size()+"; construction identities, not caller calibration");
	}
	private static float quantile(FloatList sorted, int percent){
		assert(sorted.size>0 && percent>=0 && percent<=100):"Order statistics require nonempty validated construction members";
		return sorted.array[(int)((sorted.size-1L)*percent/100)];
	}
	private float clamp(float value){return Math.max(floor, Math.min(.80f, value));}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {Member(int n){length=n;}final int length;int added=-1;String addedLabel;boolean seen;}
	private static final class Model {
		Model(Read r){ref=r;}
		void add(int length, float identity){assert(length>0):"Member extrema constrain emitted spans";values.add(identity);min=Math.min(min, length);max=Math.max(max, length);}
		final Read ref;final FloatList values=new FloatList();int min=Integer.MAX_VALUE, max;boolean checked, singleton;
	}
	private final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class);
	private final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
	private final ArrayList<Model> models=new ArrayList<Model>();
	private float floor=.60f;
	public static final String QUALITY_HEADER="locus_id\tlength\tbuilder_status\tassigned_model\tassigned_identity\tbest_model\tbest_shorter_query_identity\tat_report_cut";
	public static final String ADDITIONS_HEADER="locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction";
}
