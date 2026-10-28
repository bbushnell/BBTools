package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectIntMap;
import map.ObjectMap;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import structures.IntList;

/** Reports biological composition and lengths from builder and added memberships.
 * A member may support its original model and a separately recorded singleton.
 * These are construction memberships, not best-match counts or HBM vote depths.
 * @author Brian, Raiden
 */
public final class NcrnaModelMembershipReport {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaModelMembershipReport.class, false);
		String input=null, refs=null, assignments=null, additions=null, out=null, family=null, labels=null, prefix=null;
		int expected=-1;double purity=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){input=v;}else if(k.equals("ref")){refs=v;}else if(k.equals("assignments")){assignments=v;}
			else if(k.equals("additions")){additions=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("family")){family=v;}
			else if(k.equals("labels")){labels=v;}else if(k.equals("modelprefix")){prefix=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("purity")){purity=Double.parseDouble(v);}
			else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(input!=null && refs!=null && assignments!=null && out!=null && family!=null && labels!=null && prefix!=null
			&& !prefix.isEmpty() && expected>0 && purity>0 && purity<=1,
			"Bind complete input, model list, original assignments, emitted family, biological labels, prefix, denominator and purity threshold");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{new NcrnaModelMembershipReport().run(input, refs, assignments, additions, out, family, labels.split(",", -1), prefix, expected, purity);}
		finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Methods           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(String input, String refs, String assignments, String additions, String out, String family, String[] labels, String prefix, int expected, double purity){
		final ObjectIntMap<String> biological=new ObjectIntMap<String>(String.class), modelIds=new ObjectIntMap<String>(String.class);
		for(int i=0; i<labels.length; i++){
			require(labels[i].matches("[A-Za-z0-9_]+") && !biological.contains(labels[i]), "Biological labels must be unique TSV-safe names");biological.put(labels[i], i);
		}
		require(biological.contains(family), "The emitted family must have an explicit biological-label column");
		final ArrayList<Model> models=new ArrayList<Model>();
		RrnaResourceIO.stream(refs, r->{
			final String id=RrnaResourceIO.first(r.id);require(id.startsWith(prefix) && !modelIds.contains(id) && r.bases!=null && r.length()>0,
				"Models must have unique names under the declared emitted-family prefix: "+id);
			modelIds.put(id, models.size());models.add(new Model(r, labels.length));
		});
		require(!models.isEmpty(), "A nonempty model library is required");
		final ObjectMap<String,Member> members=assignments(assignments, family, modelIds, models, expected);
		if(additions!=null){additions(additions, family, biological, modelIds, models, members);}
		final int[] states=new int[4];final IntList groupLengths=new IntList();
		RrnaResourceIO.stream(input, r->{
			final String id=RrnaResourceIO.first(r.id);final Member m=members.get(id);
			require(m!=null && !m.seen && r.id.equals(m.header) && r.bases!=null && r.length()>0,
				"Every original assignment must retain its full input header and occur exactly once: "+id);
			final String label=RrnaResourceIO.field(r.id, "member_family");
			require(biological.contains(label) && RrnaResourceIO.field(r.id, "build_group").equals(family), "Input biological label or operational group is undeclared: "+id);
			final int bi=biological.get(label);m.seen=true;groupLengths.add(r.length());states[m.state]++;
			if(m.model>=0){models.get(m.model).add(r.length(), bi);}
			if(m.added>=0){
				final Model model=models.get(m.added);
				require(label.equals(m.addedLabel) && r.length()==m.addedLength && Arrays.equals(r.bases, model.ref.bases),
					"An added singleton must retain the exact input sequence, length and biological label: "+id);
				model.add(r.length(), bi);states[3]++;
			}
		});
		require(groupLengths.size==expected && states[0]+states[1]+states[2]==expected,
			"Builder dispositions must conserve the entire population before reporting model composition");
		long memberships=0;for(Model model:models){memberships+=model.lengths.size;}
		require(memberships==(long)states[0]+states[3], "Construction membership events equal original kept members plus explicit added contributors");
		write(out, family, labels, biological.get(family), expected, purity, models, groupLengths, states, memberships);
	}

	private static ObjectMap<String,Member> assignments(String path, String family, ObjectIntMap<String> ids, ArrayList<Model> models, int expected){
		final ObjectMap<String,Member> members=new ObjectMap<String,Member>(String.class, Member.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==6, "Builder assignment schema changed");
				final String header=p.parseString(0), id=RrnaResourceIO.first(header), label=p.parseString(4);
				final int state=p.termEquals("kept", 5)?0:p.termEquals("small", 5)?1:p.termEquals("orphan", 5)?2:-1;
				require(state>=0 && !members.contains(id) && p.termStartsWith(family+"_", 3), "Unique family-namespaced builder assignments are required: "+id);
				int model=-1;
				if(state==0){
					final String name=RrnaResourceIO.first(label);require(ids.contains(name), "Original kept assignment names a missing model: "+name);
					model=ids.get(name);models.get(model).originalMembers++;
				}else{require(label.isEmpty(), "Unassigned members must not claim a model");}
				members.put(id, new Member(header, model, state));
			}
		}finally{require(!in.close(), "Assignment read failed");}
		require(members.size()==expected, "Original assignments must include all admitted members");return members;
	}

	private static void additions(String path, String family, ObjectIntMap<String> biological, ObjectIntMap<String> ids, ArrayList<Model> models, ObjectMap<String,Member> members){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "locus_id\tmodel_id\temitted_family\tmember_family\tlength\tconstruction", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==6 && p.termEquals(family, 2) && p.termEquals("EXACT_SINGLE_MEMBER_NO_REASSIGNMENT", 5), "Added membership must preserve the declared family and singleton recipe");
				final Member m=members.get(p.parseString(0));final String name=p.parseString(1), label=p.parseString(3);
				require(m!=null && m.added<0 && ids.contains(name) && biological.contains(label) && p.parseInt(4)>0,
					"Each added singleton must name an existing input, model and biological label exactly once");
				final int index=ids.get(name);final Model model=models.get(index);
				require(!model.added && model.originalMembers==0, "An added singleton must not overwrite or double-count an original model: "+name);
				model.added=true;m.added=index;m.addedLabel=label;m.addedLength=p.parseInt(4);
			}
		}finally{require(!in.close(), "Added membership read failed");}
	}

	private static void write(String out, String family, String[] labels, int own, int expected, double cut, ArrayList<Model> models, IntList groupLengths, int[] states, long memberships){
		RrnaResourceIO.fresh(out, ".models.tsv", ".summary.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(out+".models.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int pure=0, mixed=0, empty=0;
		try{
			b.append("model_id\temitted_family\tconsensus_length\tconstruction_members\tmember_min\tmember_p10\tmember_median\tmember_p90\tmember_max\temitted_family_fraction\tpurity");
			for(String label:labels){b.tab().append(label).append("_members").tab().append(label).append("_fraction");}rows.print(b.nl());
			for(Model model:models){
				final int n=model.lengths.size;int count=0;for(int c:model.counts){count+=c;}
				require(count==n && (!model.added || n==1), "Model family counts must conserve actual assigned or singleton membership events");
				final double fraction=n==0?Double.NaN:model.counts[own]/(double)n;
				final String purity=n==0?"EMPTY":fraction>=cut?"PURE":"MIXED";
				if(n==0){empty++;}else if(fraction>=cut){pure++;}else{mixed++;}
				b.clear().append(RrnaResourceIO.first(model.ref.id)).tab().append(family).tab().append(model.ref.length()).tab().append(n);
				quantiles(b, model.lengths);b.tab();if(n==0){b.append("NA");}else{b.append(fraction, 6);}b.tab().append(purity);
				for(int c:model.counts){b.tab().append(c).tab();if(n==0){b.append("NA");}else{b.append(c/(double)n, 6);}}
				rows.print(b.nl());
			}
			summary.println("emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max");
			b.clear().append(family).tab().append(expected);for(int n:states){b.tab().append(n);}
			b.tab().append(memberships).tab().append(models.size()).tab().append(pure).tab().append(mixed).tab().append(empty).tab().append(cut, 6);
			quantiles(b, groupLengths);summary.print(b.nl());
		}finally{RrnaResourceIO.finish(rows, summary);}
		System.err.println("NCRNA_MODEL_MEMBERSHIP_REPORT_PASS members="+expected+" models="+models.size()+" membership_events="+memberships+"; construction membership, not best-match counts or HBM vote depths");
	}
	private static void quantiles(ByteBuilder b, IntList lengths){
		assert(lengths!=null):"Membership lengths must exist independently of consensus length";
		lengths.sort();for(int q:new int[]{0, 10, 50, 90, 100}){b.tab();if(lengths.size==0){b.append("NA");}else{b.append(lengths.array[(int)((lengths.size-1L)*q/100)]);}}
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Member {
		Member(String h, int m, int s){header=h;model=m;state=s;assert(s>=0 && s<3):"Builder dispositions index kept/small/orphan counters";}
		final String header;final int model, state;boolean seen;int added=-1, addedLength;String addedLabel;
	}
	private static final class Model {
		Model(Read r, int n){ref=r;counts=new int[n];assert(n>0):"Every model uses the same declared biological-label columns";}
		void add(int length, int label){assert(length>0 && label>=0 && label<counts.length):"Only validated member length and biological label may enter model statistics";lengths.add(length);counts[label]++;}
		final Read ref;final int[] counts;final IntList lengths=new IntList();int originalMembers;boolean added;
	}
}
