package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Taxonomy-stratified engineering split; whole source groups remain in one fold.
 * Exact compressed-genome bytes plus identical locus coordinates prove aliases.
 * Other source links block leakage but do not imply sequence equivalence.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSeedSplit {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaSeedSplit.class, false);
		String eligible=null, loci=null, groups=null, out=null, family=null;int expected=-1;
		long seed=Long.MIN_VALUE;double fraction=Double.NaN;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("eligible")){eligible=value;}else if(key.equals("loci")){loci=value;}else if(key.equals("groups")){groups=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("family") || key.equals("families")){require(family==null, "Declare the family list once");family=value;}
			else if(key.equals("expected")){expected=Parse.parseIntKMG(value);}else if(key.equals("seed")){seed=Long.parseLong(value);}
			else if(key.equals("holdout")){fraction=Double.parseDouble(value);}else{throw new IllegalArgumentException(key);}
		}
		require(eligible!=null && loci!=null && groups!=null && out!=null && family!=null && !family.isEmpty() && expected>0 && seed!=Long.MIN_VALUE && fraction>0 && fraction<1,
			"Bind full eligible IDs, source metadata/group ledger, family, denominator, seed and heldout fraction");
		Shared.setThreads(1);final NcrnaSeedSplit tool=new NcrnaSeedSplit();
		tool.loadSources(groups);tool.loadLoci(eligible, loci, family, expected);tool.split(fraction, seed);tool.write(out, fraction, seed);
		Shared.closeStream(pp.outstream);
	}

	private void loadSources(final String path){
		assert(sources.size()==0):"Source grouping must be fixed before splitting loci";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, "genome\tassembly\ttaxid\tfasta_sha80\tsource_group\tpilot_exposed_group", path);
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==6, "Source-group schema changed");
				final String name=p.parseString(0), hash=p.parseString(3), key=p.parseString(4);
				require(!sources.contains(name) && !key.isEmpty(), "Source names must be unique and group keys nonempty");
				require(hash.equals("NA") || hash.matches("[0-9a-f]{20}"), "Source content identity must use sha80 or explicit NA");
				require(p.termEquals("true", 5) || p.termEquals("false", 5), "Pilot exposure must be explicit");
				Group group=byGroup.get(key);
				if(group==null){group=new Group(key, p.termEquals("true", 5));byGroup.put(key, group);}
				require(group.pilot==p.termEquals("true", 5), "Exposure must be consistent across a transitive source group");
				sources.put(name, new Source(name, p.parseString(1), p.parseInt(2), hash, group));
			}
		}finally{require(!in.close(), "Source-group ledger read failed");}
	}

	private void loadLoci(final String eligible, final String path, final String family, final int expected){
		final ObjectSet<String> wanted=new ObjectSet<String>(String.class), seen=new ObjectSet<String>(String.class);
		final ObjectSet<String> families=new ObjectSet<String>(String.class);
		final String[] declaredFamilies=family.split(",");
		for(String name:declaredFamilies){require(!name.isEmpty() && families.add(name), "Family names must be nonempty and distinct");}
		ByteFile in=RrnaResourceIO.open(eligible);
		try{for(byte[] row;(row=in.nextLine())!=null;){
			final String id=new String(row, java.nio.charset.StandardCharsets.US_ASCII);
			boolean declared=false;for(String name:declaredFamilies){declared|=id.startsWith(name+"_");}
			require(declared && wanted.add(id), "Eligible IDs must be unique members of the declared families: "+id);
		}}finally{require(!in.close(), "Eligible-ID read failed");}
		require(wanted.size()==expected, "Full eligible population changed: "+wanted.size()+" expected="+expected);
		in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			final byte[] header=in.nextLine();require(header!=null, "Seed metadata header is required");
			final String schema=new String(header, java.nio.charset.StandardCharsets.US_ASCII);
			final boolean compact=schema.equals(NcrnaClanSeedPrep.SEED_HEADER);
			require(compact || families.size()==1 && schema.equals(NcrnaCmLocusAudit.HEADER), "Multiple families require explicit source/member/build labels in compact seed metadata");
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==(compact?12:21), "Source-locus schema changed");final String id=p.parseString(0);
				if(!wanted.contains(id)){continue;}
				require(seen.add(id), "Duplicate eligible locus in the source ledger: "+id);
				final Source source=sources.get(p.parseString(1));require(source!=null, "Eligible locus source has no fixed group: "+id);
				require(source.taxid==p.parseInt(3) && source.assembly.equals(p.parseString(2)) && !source.hash.equals("NA"), "Eligible source identity/hash disagrees with its group ledger: "+id);
				final String key=p.parseString(4)+"/"+p.parseString(5);Stratum s=byStratum.get(key);
				if(s==null){s=new Stratum(key);byStratum.put(key, s);strata.add(s);}
				final Group g=source.group;
				if(g.loci==0){groups.add(g);}
				if(g.phylumNames.add(key)){g.phyla.add(s);s.groups.add(g);}
				g.loci++;s.loci++;
				String alias=source.hash+'\t'+p.parseString(6)+'\t'+p.parseString(7)+'\t'+p.parseString(8);
				if(compact){
					final String original=p.parseString(9), member=p.parseString(10), build=p.parseString(11);
					require(families.contains(original) && id.startsWith(original+"_") && families.contains(build)
						&& (member.equals(original) && build.equals(original) || member.equals("B12clade")), "Compact source/biological/build labels must retain admitted provenance: "+id);
					// Coincident cores with distinct labels remain separate contributors.
					alias+='\t'+original+'\t'+member+'\t'+build;
				}
				loci.add(new Locus(id, source, s, alias));
			}
		}finally{require(!in.close(), "Source-locus read failed");}
		require(loci.size()==expected, "Every full-population eligible ID must have exactly one source row");
		Collections.sort(loci, (a,b)->a.id.compareTo(b.id));Collections.sort(groups, (a,b)->a.key.compareTo(b.key));
		Collections.sort(strata, (a,b)->a.key.compareTo(b.key));
		final ObjectMap<String,Locus> representatives=new ObjectMap<String,Locus>(String.class, Locus.class);
		for(Locus l:loci){
			Locus rep=representatives.get(l.alias);
			if(rep==null){rep=l;representatives.put(l.alias, l);l.source.group.representatives++;unique++;}
			require(rep.source.group==l.source.group, "Identical genome bytes must already be linked in the source-group ledger");
			l.representative=rep;
		}
	}

	private void split(final double fraction, final long seed){
		assert(!groups.isEmpty()):"Seed feasibility requires a nonempty construction population";
		final Random random=new Random(seed);
		for(Stratum s:strata){
			Collections.sort(s.groups, (a,b)->a.key.compareTo(b.key));
			final int n=s.groups.size();s.target=n<2?0:Math.max(1, Math.min(n-1, (int)Math.round(n*fraction)));
			final ArrayList<Group> candidates=new ArrayList<Group>();
			for(Group g:s.groups){if(g.phyla.size()==1){candidates.add(g);}}
			Collections.shuffle(candidates, random);
			for(int i=0;i<Math.min(s.target, candidates.size());i++){candidates.get(i).heldout=true;s.heldoutGroups++;}
			require(s.heldoutGroups<n, "Every represented phylum must retain a training source group: "+s.key);
		}
		for(Locus l:loci){if(l.source.group.heldout){l.stratum.heldoutLoci++;}}
	}

	private void write(final String prefix, final double fraction, final long seed){
		RrnaResourceIO.fresh(prefix, ".train.ids.txt", ".train.raw.ids.txt", ".holdout.ids.txt", ".loci.tsv", ".groups.tsv", ".phyla.tsv", ".summary.tsv");
		final ByteStreamWriter train=RrnaResourceIO.writer(prefix+".train.ids.txt"), holdout=RrnaResourceIO.writer(prefix+".holdout.ids.txt");
		final ByteStreamWriter trainAll=RrnaResourceIO.writer(prefix+".train.raw.ids.txt");
		final ByteStreamWriter rows=RrnaResourceIO.writer(prefix+".loci.tsv"), gg=RrnaResourceIO.writer(prefix+".groups.tsv"), pp=RrnaResourceIO.writer(prefix+".phyla.tsv"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();int trainRaw=0, heldRaw=0, trainUnique=0, heldUnique=0, heldGroups=0;
		try{
			rows.println("locus_id\tgenome\tsource_group\tstratum\trepresentative_id\tfold\tpilot_exposed_group");
			for(Locus l:loci){
				final boolean held=l.source.group.heldout, representative=l==l.representative;
				if(held){heldRaw++;holdout.println(l.id);if(representative){heldUnique++;}}
				else{trainRaw++;trainAll.println(l.id);if(representative){trainUnique++;train.println(l.id);}}
				rows.print(b.clear().append(l.id).tab().append(l.source.name).tab().append(l.source.group.key).tab().append(l.stratum.key).tab().append(l.representative.id)
					.tab().append(held?"HELD_OUT":"TRAIN").tab().append(l.source.group.pilot).nl());
			}
			require(trainRaw+heldRaw==loci.size() && trainUnique+heldUnique==unique && trainUnique>0 && heldRaw>0,
				"Split must conserve raw/alias-collapsed counts and provide nonempty construction/evaluation folds");
			gg.println("source_group\tfold\tloci\trepresentatives\tphyla\tpilot_exposed_group\treason");
			for(Group g:groups){
				if(g.heldout){heldGroups++;}
				final String reason=g.heldout?"STRATIFIED_HOLDOUT":g.phyla.size()>1?"MULTIPHYLUM_GROUP_TRAIN":g.phyla.get(0).groups.size()==1?"SINGLETON_PHYLUM_TRAIN":"STRATIFIED_TRAIN";
				gg.print(b.clear().append(g.key).tab().append(g.heldout?"HELD_OUT":"TRAIN").tab().append(g.loci).tab().append(g.representatives).tab().append(g.phyla.size())
					.tab().append(g.pilot).tab().append(reason).nl());
			}
			pp.println("stratum\tsource_groups\ttrain_groups\theldout_groups\ttarget_heldout_groups\tloci\ttrain_loci\theldout_loci\texception");
			for(Stratum s:strata){
				final String reason=s.groups.size()==1?"SINGLETON_TRAIN_ONLY":s.heldoutGroups<s.target?"MULTIPHYLUM_GROUP_CONSTRAINT":"NONE";
				pp.print(b.clear().append(s.key).tab().append(s.groups.size()).tab().append(s.groups.size()-s.heldoutGroups).tab().append(s.heldoutGroups).tab().append(s.target)
					.tab().append(s.loci).tab().append(s.loci-s.heldoutLoci).tab().append(s.heldoutLoci).tab().append(reason).nl());
			}
			summary.println("eligible_loci\tunique_source_loci\talias_loci\ttrain_loci_raw\ttrain_representatives\theldout_loci_raw\theldout_representatives\tsource_groups\ttrain_groups\theldout_groups\ttarget_group_fraction\tseed\tinterpretation");
			summary.print(b.clear().append(loci.size()).tab().append(unique).tab().append(loci.size()-unique).tab().append(trainRaw).tab().append(trainUnique).tab().append(heldRaw).tab().append(heldUnique)
				.tab().append(groups.size()).tab().append(groups.size()-heldGroups).tab().append(heldGroups).tab().append(fraction, 4).tab().append(seed).append("\tTAXONOMY_STRATIFIED_SEED_FEASIBILITY\n"));
		}finally{RrnaResourceIO.finish(train, trainAll, holdout, rows, gg, pp, summary);}
		System.err.println("NCRNA_SEED_SPLIT_PASS eligible="+loci.size()+" trainRepresentatives="+trainUnique+" heldoutRaw="+heldRaw+"; no sequence-cluster protection or biological-completeness claim");
	}

	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Source {
		Source(String n, String a, int t, String h, Group g){name=n;assembly=a;taxid=t;hash=h;group=g;}
		final String name, assembly, hash;final int taxid;final Group group;
	}
	private static final class Group {
		Group(String k, boolean p){key=k;pilot=p;}final String key;final boolean pilot;boolean heldout;int loci, representatives;
		final ObjectSet<String> phylumNames=new ObjectSet<String>(String.class);final ArrayList<Stratum> phyla=new ArrayList<Stratum>();
	}
	private static final class Stratum {
		Stratum(String k){key=k;}final String key;final ArrayList<Group> groups=new ArrayList<Group>();int loci, heldoutLoci, heldoutGroups, target;
	}
	private static final class Locus {
		Locus(String i, Source s, Stratum p, String a){id=i;source=s;stratum=p;alias=a;}final String id, alias;final Source source;final Stratum stratum;Locus representative;
	}
	private final ObjectMap<String,Source> sources=new ObjectMap<String,Source>(String.class, Source.class);
	private final ObjectMap<String,Group> byGroup=new ObjectMap<String,Group>(String.class, Group.class);
	private final ObjectMap<String,Stratum> byStratum=new ObjectMap<String,Stratum>(String.class, Stratum.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<Group> groups=new ArrayList<Group>();private final ArrayList<Stratum> strata=new ArrayList<Stratum>();private int unique;
}
