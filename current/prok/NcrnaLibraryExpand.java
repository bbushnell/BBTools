package prok;

import java.nio.charset.StandardCharsets;
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
import stream.Read;
import structures.ByteBuilder;

/** Adds exact singleton models for uncovered admitted members; preserves every old model.
 * The original builder assignments stay immutable. Added-model membership is separate.
 * Optional inputfamily= declares one producing-CM family for native source IDs;
 * it must equal family= and every ID must use that family's input namespace.
 * Existing biological/group tags must agree. Without it, historical tag checks apply.
 * @author Brian, Raiden
 */
public final class NcrnaLibraryExpand {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaLibraryExpand.class, false);
		String in=null, ref=null, hbm=null, report=null, out=null, family=null, prefix=null;
		String inputFamily=null;
		int expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');
			require(eq>0, "Use flag=value");
			final String a=arg.substring(0, eq).toLowerCase(Locale.ROOT), b=arg.substring(eq+1);
			if(a.equals("in")){in=b;}
			else if(a.equals("ref")){ref=b;}
			else if(a.equals("hbm")){hbm=b;}
			else if(a.equals("report")){report=b;}
			else if(a.equals("outprefix")){out=b;}
			else if(a.equals("family")){family=b;}
			else if(a.equals("prefix")){prefix=b;}
			else if(a.equals("expected")){expected=Integer.parseInt(b);}
			else if(a.equals("inputfamily")){inputFamily=b;}
			else{throw new IllegalArgumentException("Unknown option: "+a);}
		}
		require(in!=null && ref!=null && hbm!=null && report!=null && out!=null && family!=null && prefix!=null && !prefix.isEmpty() && expected>0,
			"Bind the entire admitted corpus, prior library and exhaustive report, output, group, new model prefix and population");
		require(inputFamily==null || inputFamily.equals(family),
			"Explicit input family must equal the emitted family; cross-family routing requires original metadata");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{new NcrnaLibraryExpand().run(in, ref, hbm, report, out, family, prefix, expected, inputFamily);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(String input, String refs, String hbm, String report, String out, String family,
		String prefix, int expected, String inputFamily){
		final ArrayList<Read> library=new ArrayList<Read>();RrnaResourceIO.stream(refs, library::add);
		final ArrayList<NucleotideHbmCodec.Model> models=new ArrayList<NucleotideHbmCodec.Model>();
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		try(NucleotideHbmCodec.Reader reader=new NucleotideHbmCodec.Reader(hbm)){
			for(NucleotideHbmCodec.Model model; (model=reader.next())!=null; ){
				require(models.size()<library.size(), "HBM contains more models than consensus FASTA");
				final Read r=library.get(models.size());final String name=RrnaResourceIO.first(r.id);
				require(name.equals(RrnaResourceIO.first(new String(model.header, 1, model.header.length-1, StandardCharsets.UTF_8)))
					&& Arrays.equals(r.bases, model.bases) && names.add(name), "Prior FASTA/HBM names, order and bases must agree: "+name);
				models.add(model);
			}
		}
		require(!models.isEmpty() && models.size()==library.size(), "A nonempty concordant original library is required");
		final ObjectMap<String,Member> members=loadReport(report, names, expected);
		final ArrayList<Read> additions=new ArrayList<Read>();final int[] count={0};
		RrnaResourceIO.stream(input, r->{
			final String id=RrnaResourceIO.first(r.id);final Member m=members.get(id);
			require(m!=null && !m.seen && r.bases!=null && r.length()==m.length && r.length()>0,
				"Every reported member must occur once in the complete unchanged-length admitted FASTA: "+id);
			validateInputFamily(r.id, family, inputFamily);
			m.seen=true;count[0]++;
			if(!m.covered){
				for(byte base:r.bases){require(RrnaResourceIO.base(base)>=0, "Singleton models require unambiguous admitted DNA: "+id);}
				if(inputFamily==null){RrnaResourceIO.field(r.id, "member_family");}
				require(names.add(prefix+id), "Expanded model name collides with an existing model: "+id);additions.add(r);
			}
		});
		require(count[0]==expected, "Complete input/report membership must be conserved before expansion: "+count[0]+" expected="+expected);
		RrnaResourceIO.fresh(out, ".consensus.fa", ".models.hbmc", ".additions.tsv");
		final ByteStreamWriter fasta=RrnaResourceIO.writer(out+".consensus.fa"), ledger=RrnaResourceIO.writer(out+".additions.tsv");
		final ByteBuilder b=new ByteBuilder();
		try(NucleotideHbmCodec.Writer writer=new NucleotideHbmCodec.Writer(out+".models.hbmc", true, false)){
			for(int i=0; i<models.size(); i++){
				writer.add(models.get(i));final Read r=library.get(i);fasta.print(b.clear().append('>').append(r.id).nl().append(r.bases).nl());
			}
			ledger.println("locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction");
			for(Read r:additions){
				final String id=RrnaResourceIO.first(r.id), model=prefix+id, label=model+" n=1";
				final int[] counts=new int[Math.multiplyExact(r.length(), 5)];
				for(int i=0; i<r.length(); i++){counts[5*i+RrnaResourceIO.base(r.bases[i])]=1;}
				// Same singleton counts as TrnaConsensusBuilder.buildConsensusAndGraph:
				// BaseGraph's constructor gives the exact reference one vote per position.
				writer.add(new NucleotideHbmCodec.Model((">"+label+"\tdepth=1").getBytes(StandardCharsets.UTF_8), r.bases, counts));
				fasta.print(b.clear().append('>').append(label).nl().append(r.bases).nl());
				final String biological=inputFamily==null ? RrnaResourceIO.field(r.id, "member_family") : inputFamily;
				ledger.print(b.clear().append(id).tab().append(model).tab().append(family).tab().append(biological)
					.tab().append(r.length()).append("\tEXACT_SINGLE_MEMBER_NO_REASSIGNMENT\n"));
			}
			writer.finish();
		}finally{RrnaResourceIO.finish(fasta, ledger);}
		System.err.println("NCRNA_LIBRARY_EXPAND_PASS original_models="+models.size()+" added_models="+additions.size()+" members="+expected+
			"; original assignments unchanged, added memberships separate; exhaustive remeasurement required");
	}

	/** Validates the historical tag contract or the explicit single-producing-CM binding.
	 * Never changes a header, sequence, member disposition, or model contribution.
	 */
	private static void validateInputFamily(String header, String family, String inputFamily){
		assert(header!=null && family!=null):"Expansion requires intact source headers and an emitted group";
		if(inputFamily==null){
			require(RrnaResourceIO.field(header, "build_group").equals(family),
				"Foreign operational group in admitted input: "+RrnaResourceIO.first(header));
			return;
		}
		require(RrnaResourceIO.first(header).startsWith(inputFamily+"_input_"),
			"Explicit input family must match every original source namespace: "+header);
		final String group=RrnaResourceIO.optionalField(header, "build_group");
		final String biological=RrnaResourceIO.optionalField(header, "member_family");
		require((group==null || group.equals(inputFamily)) && (biological==null || biological.equals(inputFamily)),
			"Existing metadata conflicts with the explicit producing-CM family: "+header);
	}

	private static ObjectMap<String,Member> loadReport(String path, ObjectSet<String> names, int expected){
		final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "locus_id\tlength\tbuilder_status\tassigned_model\tassigned_identity\tbest_model\tbest_shorter_query_identity\tat_report_cut", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==8, "Native exhaustive member report schema changed");
				final String id=p.parseString(0);final int length=p.parseInt(1);final double best=p.parseDouble(6);
				require(!members.contains(id) && length>0 && names.contains(p.parseString(5)) && Double.isFinite(best) && best>=0 && best<=1
					&& (p.termEquals("true", 7) || p.termEquals("false", 7)), "Unique valid member rows and known best models are required: "+id);
				members.put(id, new Member(length, p.termEquals("true", 7)));
			}
		}finally{require(!in.close(), "Exhaustive member report read failed");}
		require(members.size()==expected, "Every admitted member, including the covered members, must be in the report");return members;
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(int n, boolean c){length=n;covered=c;assert(n>0):"The report must retain a real positive-length core for singleton construction";}
		final int length;final boolean covered;boolean seen;
	}
}
