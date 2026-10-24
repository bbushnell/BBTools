package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.ReadWrite;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Descriptive locus construction, never truth admission.
 * G11's 2026-10-18 Q3 rule: partition by assembly/native seqid/strand/family;
 * seed class GA+full, GA+partial, below-GA precedes bits DESC, genomic length
 * DESC, accession/start/end/run/table/line ASC (G11 revision02:11UTC).
 * Subtype is the highest-bit GA member, independently of the geometry seed.
 * A spanning hit positively overlaps two GA-passing, full-model neighbors whose
 * mutual overlap is less than half the shorter (G11 revision01:19UTC).
 * Partial/below-GA pieces cannot demote a containing hit. Non-spanning seeds
 * and members are processed first.
 * Deferred spanning hits join the earliest qualifying seed, or seed a residual.
 * Membership requires at least half the MEMBER length; seed spans never grow.
 * Every original overlay row survives in the members or invalid companion.
 * @author Raiden */
public final class EukPanelLoci {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, EukPanelLoci.class, false);
		String in=null, models=null, scans=null, out=null;long expected=-1;int expectedModels=-1, expectedScans=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Expected flag=value");
			final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("models")){models=v;}else if(k.equals("scans")){scans=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedhits")){expected=Long.parseLong(v);}
			else if(k.equals("expectedmodels")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("expectedscans")){expectedScans=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(in!=null && models!=null && scans!=null && out!=null && expected>=0 && expectedModels>0 && expectedScans>0, "Explicit inputs and denominators required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelLoci x=new EukPanelLoci();x.loadModels(models, expectedModels);x.loadScans(scans, expectedScans);
		x.load(in, expected);x.construct();x.write(out);Shared.closeStream(pp.outstream);
	}
	void loadModels(String path, int expected){
		assert(models.isEmpty()):"One pinned model roster defines fold membership";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, "bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80", path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==8, "Model schema changed");
				final Model m=new Model(p.parseString(0), p.parseString(1), p.parseString(2), p.parseInt(4), p.parseDouble(5));
				require(SeedOffsetTable.validSha80(m.bundle) && m.clen>0 && Double.isFinite(m.ga), "Invalid model metadata");
				require(models.put(m.key(), m)==null, "Duplicate model identity");
			}
		}finally{require(!in.close(), "Model read failed");}
		require(models.size()==expected, "Model denominator changed");
	}
	void loadScans(String path, int expected){
		assert(scans.isEmpty()):"Zero-hit scans supply reporting denominators";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, "selection_bin\t"+EukPanelCmProposal.SUMMARY_HEADER, path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==14, "Scan schema changed");
				final String bin=p.parseString(0), assembly=p.parseString(1), run=p.parseString(2), accession=p.parseString(3), bundle=p.parseString(4);
				final Model m=models.get(bundle+'\t'+accession);require(m!=null, "Scan has unknown model");
				final String key=assembly+'\t'+run+'\t'+bundle+'\t'+accession;
				final Scan s=new Scan(bin, p.parseLong(5));require(s.expected>=0 && scans.put(key, s)==null, "Duplicate or invalid scan");
				final String group=bin+'\t'+m.family;if(counts.get(group)==null){counts.put(group, new long[6]);}
			}
		}finally{require(!in.close(), "Scan read failed");}
		require(scans.size()==expected, "Scan denominator changed");
	}
	void load(String path, long expected){
		assert(hits.isEmpty()):"Every raw hit is loaded once before any assignments";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in, HEADER, path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==40 && p.termEquals("REVIEW_ONLY",39), "Only the accepted descriptive overlay is input");
			final Hit h=new Hit(line, p);final Model m=models.get(h.h.bundle+'\t'+h.h.model);
			require(m!=null && m.family.equals(h.h.family) && m.clen==h.h.clen && m.ga==p.parseDouble(26), "Hit model metadata disagrees with pinned roster");
			final Scan s=scans.get(h.h.summaryKey);require(s!=null && s.bin.equals(h.h.bin), "Hit has no matching bin/scan denominator");s.observed++;
			require(seen.add(h.h.key), "Duplicate run/line identity");
			require(p.termEquals(h.h.full?"true":"false",32) && p.termEquals(h.h.gaPass?"true":"false",35), "Overlay span/GA flag disagrees with raw metrics");
			if(h.h.valid){
				require(!h.h.nativeSeqid.equals("NA") && h.h.contigLength>=h.h.stop, "Valid interval must resolve within its contig");
				final String edge=h.h.start==1?(h.h.stop==h.h.contigLength?"BOTH":"LEFT"):h.h.stop==h.h.contigLength?"RIGHT":"INTERNAL";
				require(h.h.edge.equals(edge), "Contig contact disagrees with physical coordinates");
				ArrayList<Hit> group=partitions.get(h.partition);if(group==null){group=new ArrayList<Hit>();partitions.put(h.partition, group);}group.add(h);
			}else{require(h.h.edge.equals("UNASSESSABLE"), "Invalid interval cannot acquire physical edge evidence");invalid++;}
			hits.add(h);
		}}finally{require(!in.close(), "Overlay read failed");}
		require(hits.size()==expected, "Raw-hit denominator changed: "+hits.size());
		for(Scan s:scans.values()){if(s!=null){require(s.observed==s.expected, "Hit count differs from its scan denominator");}}
	}
	void construct(){
		assert(loci.isEmpty()):"Locus IDs are assigned once from sorted partitions";
		final ArrayList<String> keys=new ArrayList<String>();for(String k:partitions.keys()){if(k!=null){keys.add(k);}}Collections.sort(keys);
		for(String key:keys){buildPartition(partitions.get(key), loci);}
		long n=0;for(Locus locus:loci){
			require(locus.seed.locus==locus && !locus.members.isEmpty(), "Every locus must own its seed");
			for(Hit h:locus.members){require(h.locus==locus && h.partition.equals(locus.seed.partition) && qualifies(h, locus.seed), "Membership crossed identity or seed-overlap rule");n++;}
			findAmbiguity(locus);
		}
		require(n+invalid==hits.size(), "Locus members plus invalid companions must conserve every input hit");
	}
	public static void buildPartition(ArrayList<Hit> group, ArrayList<Locus> out){
		assert(!group.isEmpty()):"Only populated partitions enter the overlap sweep";
		for(Hit h:group){require(h.locus==null && h.neighbors.isEmpty(), "Partition may be constructed only once");}
		Collections.sort(group, POSITION);
		for(int i=0; i<group.size(); i++){
			final Hit a=group.get(i);
			for(int j=i+1; j<group.size(); j++){
				final Hit b=group.get(j);if(b.h.start>a.h.stop){break;}a.neighbors.add(b);b.neighbors.add(a);
			}
		}
		for(Hit h:group){
			for(int i=0; i<h.neighbors.size() && !h.spanning; i++){
				for(int j=i+1; j<h.neighbors.size(); j++){
					final Hit a=h.neighbors.get(i), b=h.neighbors.get(j);
					if(substantial(a) && substantial(b) && separate(a,b)){h.spanning=true;h.witnessA=a;h.witnessB=b;break;}
				}
			}
		}
		Collections.sort(group, RANK);
		for(Hit h:group){if(!h.spanning && h.locus==null){
			final Locus locus=new Locus(h, out.size()+1);out.add(locus);
			for(Hit member:h.neighbors){if(!member.spanning && member.locus==null && qualifies(member,h)){locus.add(member);}}
		}}
		for(Hit h:group){if(h.spanning){
			Locus best=null;
			for(Hit neighbor:h.neighbors){
				final Locus candidate=neighbor.locus;
				if(candidate!=null && candidate.seed==neighbor && qualifies(h,neighbor) && (best==null || candidate.id<best.id)){best=candidate;}
			}
			if(best==null){best=new Locus(h,out.size()+1);out.add(best);}else{best.add(h);}
		}}
		for(Hit h:group){require(h.locus!=null, "Every valid hit needs exactly one locus");}
	}
	public static void findAmbiguity(Locus locus){
		assert(locus.seed!=null):"Diagnostic flags refer to a fixed seed and its assigned members";
		for(int i=0; i<locus.members.size() && locus.witnessA==null; i++){
			final Hit a=locus.members.get(i);if(a==locus.seed){continue;}
			for(int j=i+1; j<locus.members.size(); j++){
				final Hit b=locus.members.get(j);if(b!=locus.seed && substantial(a) && substantial(b) && separate(a,b)){locus.witnessA=a;locus.witnessB=b;break;}
			}
		}
	}
	void write(String prefix){
		assert(!counts.isEmpty()):"Output includes all reporting bins, even folds with zero hits";
		RrnaResourceIO.fresh(prefix, ".loci.tsv.gz", ".locus_members.tsv.gz", ".invalid_hits.tsv.gz", ".locus_counts.tsv", ".member_distribution.tsv", ".locus_cases.tsv");
		final ByteStreamWriter rows=RrnaResourceIO.writer(prefix+".loci.tsv.gz"), members=RrnaResourceIO.writer(prefix+".locus_members.tsv.gz"), bad=RrnaResourceIO.writer(prefix+".invalid_hits.tsv.gz");
		final ByteStreamWriter totals=RrnaResourceIO.writer(prefix+".locus_counts.tsv"), dist=RrnaResourceIO.writer(prefix+".member_distribution.tsv"), cases=RrnaResourceIO.writer(prefix+".locus_cases.tsv");
		final ObjectMap<String,long[]> distribution=new ObjectMap<String,long[]>(String.class,long[].class);final ByteBuilder b=new ByteBuilder();
		long assigned=0, spanning=0, ambiguous=0;
		try{
			rows.println(LOCUS_HEADER);
			members.println(HEADER+"\tlocus_id\trole\tspanning_hit");bad.println(HEADER);
			cases.println("locus_id\treason\tassembly\tnative_seqid\tstrand\tfold_family\tseed_run\tseed_line\tseed_start\tseed_end\twitness_a_run\twitness_a_line\twitness_a_start\twitness_a_end\twitness_b_run\twitness_b_line\twitness_b_start\twitness_b_end");
			for(Hit h:hits){if(!h.h.valid){bad.println(h.h.raw);continue;}
				members.print(b.clear().append(h.h.raw).tab().append(h.locus.id).tab().append(h==h.locus.seed?"SEED":"MEMBER").tab().append(h.spanning?1:0).nl());assigned++;
			}
			for(Locus locus:loci){
				final Hit seed=locus.seed;final EukPanelCmProposal.Hit h=seed.h;
				int ns=0;boolean ga=false;final ObjectSet<String> accessions=new ObjectSet<String>(String.class);
				final ArrayList<String> names=new ArrayList<String>();
				for(Hit member:locus.members){ga|=member.h.gaPass;if(member.spanning){ns++;}if(accessions.add(member.h.model)){names.add(member.h.model);}}
				Collections.sort(names);final boolean flag=seed.spanning || locus.witnessA!=null;spanning+=ns;if(flag){ambiguous++;}
				b.clear().append(locus.id).tab().append(h.assembly).tab().append(h.nativeSeqid).tab().append(h.strand).tab().append(h.start).tab().append(h.stop).tab().append(h.family)
					.tab().append(h.bin).tab().append(h.compartment).tab().append(h.model).tab().append(h.run).tab().append(seed.table).tab().append(h.lineNumber).tab().append(seed.bitsText)
					.tab().append(ga?1:0).tab().append(h.full?1:0).tab().append(h.trunc.equals("no")?0:1).tab().append(h.edge).tab().append(locus.members.size()).tab();
				for(int i=0; i<names.size(); i++){if(i>0){b.append(';');}b.append(names.get(i));}
				final Hit subtype=subtype(locus.members);
				rows.print(b.tab().append(seed.spanning?1:0).tab().append(ns).tab().append(flag?"SPANNING_AMBIGUOUS":"NONE").tab().append("DESCRIPTIVE_ONLY").tab().append(h.gaPass?1:0).tab().append(subtype==null?"NA":subtype.h.model).nl());
				final String key=h.bin+'\t'+h.family;final long[] c=counts.get(key);require(c!=null, "Unknown reporting bin/fold");c[0]++;c[1]+=locus.members.size();if(ga){c[2]++;}else{c[3]++;}c[4]+=ns;if(flag){c[5]++;}
				EukPanelCmProposal.increment(distribution, key+'\t'+locus.members.size(), 1);
				if(seed.spanning){writeCase(cases,b,locus,"SPANNING_SEED",seed.witnessA,seed.witnessB);}
				if(locus.witnessA!=null){writeCase(cases,b,locus,"SEPARATE_MEMBERS",locus.witnessA,locus.witnessB);}
			}
			totals.println("selection_bin\tfold_family\tloci\tmembers\tga_supported_loci\tbelow_ga_loci\tspanning_hits\tambiguous_loci");
			long total=0;for(String k:EukPanelCmProposal.sorted(counts)){b.clear().append(k);for(long n:counts.get(k)){b.tab().append(n);}totals.print(b.nl());total+=counts.get(k)[1];}
			require(total==assigned && assigned+invalid==hits.size(), "Published counts must conserve all original hits");
			dist.println("selection_bin\tfold_family\tmembers_per_locus\tloci");for(String k:EukPanelCmProposal.sorted(distribution)){dist.print(b.clear().append(k).tab().append(distribution.get(k)[0]).nl());}
		}finally{RrnaResourceIO.finish(rows,members,bad,totals,dist,cases);}
		System.err.println("LOCUS_DIAGNOSTIC_COMPLETE hits="+hits.size()+" members="+assigned+" invalid="+invalid+" loci="+loci.size()+" spanning_hits="+spanning+" ambiguous_loci="+ambiguous+" no_truth_admission=true");
	}
	static void writeCase(ByteStreamWriter out, ByteBuilder b, Locus locus, String reason, Hit a, Hit z){
		assert(a!=null && z!=null):"Every geometry warning has two retained witness hits";
		final EukPanelCmProposal.Hit h=locus.seed.h;
		b.clear().append(locus.id).tab().append(reason).tab().append(h.assembly).tab().append(h.nativeSeqid).tab().append(h.strand).tab().append(h.family)
			.tab().append(h.run).tab().append(h.lineNumber).tab().append(h.start).tab().append(h.stop);
		for(Hit x:new Hit[]{a,z}){b.tab().append(x.h.run).tab().append(x.h.lineNumber).tab().append(x.h.start).tab().append(x.h.stop);}out.print(b.nl());
	}
	public static boolean qualifies(Hit member, Hit seed){assert(member.partition.equals(seed.partition)):"Overlap alone must never merge distinct families or strands";return overlap(member,seed)>=((member.length()+1)/2);}
	static boolean substantial(Hit h){return h.h.gaPass && h.h.full;}
	static int rankClass(Hit h){return h.h.gaPass?(h.h.full?0:1):2;}
	public static Hit subtype(ArrayList<Hit> members){
		assert(!members.isEmpty()):"Subtype evidence comes from the conserved members of one locus";
		Hit best=null;for(Hit h:members){if(h.h.gaPass && (best==null || BIT_RANK.compare(h,best)<0)){best=h;}}return best;
	}
	static boolean separate(Hit a, Hit b){assert(a.partition.equals(b.partition)):"Spanning witnesses must belong to the same biological fold partition";return overlap(a,b)<((Math.min(a.length(),b.length())+1)/2);}
	static long overlap(Hit a, Hit b){return Math.max(0, Math.min(a.h.stop,b.h.stop)-Math.max(a.h.start,b.h.start)+1);}
	public static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final class Hit{
		public Hit(byte[] raw, LineParser1 p){this(raw,p,null);}
		/** Explicit SRP union changes grouping only; original model/family/raw bytes survive. */
		public Hit(byte[] raw, LineParser1 p,String union){h=new EukPanelCmProposal.Hit(raw,p);h.bin=p.parseString(29);h.contigLength=p.termEquals("NA",30)?-1:p.parseLong(30);h.edge=p.parseString(31);
			require(union==null || union.equals("SRP_UNION") && (h.family.equals("SRP_small") || h.family.equals("SRP_large")),"Only the reviewed SRP physical fold may override original grouping");
			bits=p.parseDouble(24);bitsText=p.parseString(24);table=p.parseString(2);partition=h.assembly+'\t'+h.nativeSeqid+'\t'+h.strand+'\t'+(union==null?h.family:union);}
		long length(){return h.stop-h.start+1;}
		public final EukPanelCmProposal.Hit h;final String table,partition,bitsText;final double bits;
		final ArrayList<Hit> neighbors=new ArrayList<Hit>();public boolean spanning;Hit witnessA,witnessB;public Locus locus;
	}
	public static final class Locus{
		public Locus(Hit h, int id_){seed=h;id=id_;add(h);}
		public void add(Hit h){require(h.locus==null, "A hit cannot belong to two loci");h.locus=this;members.add(h);}
		public final Hit seed;final int id;public final ArrayList<Hit> members=new ArrayList<Hit>();public Hit witnessA,witnessB;
	}
	static final class Model{
		Model(String b,String a,String f,int c,double g){bundle=b;accession=a;family=f;clen=c;ga=g;}
		String key(){return bundle+'\t'+accession;}
		final String bundle,accession,family;final int clen;final double ga;
	}
	static final class Scan{Scan(String b,long n){bin=b;expected=n;}final String bin;final long expected;long observed;}
	static final Comparator<Hit> BIT_RANK=(a,b)->{
		int c=Double.compare(b.bits,a.bits);if(c==0){c=Long.compare(b.length(),a.length());}if(c==0){c=a.h.model.compareTo(b.h.model);}
		if(c==0){c=Long.compare(a.h.start,b.h.start);}if(c==0){c=Long.compare(a.h.stop,b.h.stop);}if(c==0){c=a.h.run.compareTo(b.h.run);}
		if(c==0){c=a.table.compareTo(b.table);}if(c==0){c=Long.compare(a.h.lineNumber,b.h.lineNumber);}return c;
	};
	public static final Comparator<Hit> RANK=(a,b)->{final int c=Integer.compare(rankClass(a),rankClass(b));return c==0?BIT_RANK.compare(a,b):c;};
	static final Comparator<Hit> POSITION=(a,b)->{final int c=Long.compare(a.h.start,b.h.start);return c==0?RANK.compare(a,b):c;};
	public static final String LOCUS_HEADER="locus_id\tassembly\tnative_seqid\tstrand\tstart\tend\tfold_family\tselection_bin\tactual_compartment\tseed_accession\tseed_run\tseed_table\tseed_line\tseed_bits\tga_pass_any\tfull_model_span\tinfernal_truncated\tcontig_edge_contact\tmembers\tmember_accessions\tspanning_seed\tspanning_members\tflags\tstatus\tseed_ga_pass\tsubtype_accession";
	public static final String HEADER=EukPanelCmProposal.HIT_HEADER+"\tselection_bin\tcontig_length\tcontig_edge_contact\tfull_model_span\tmissing_model_5prime\tmissing_model_3prime\tga_pass\tsame_model_overlaps\tsame_family_other_model_overlaps\tcross_family_overlaps\tproposal_status";
	final ObjectMap<String,Model> models=new ObjectMap<String,Model>(String.class,Model.class);
	final ObjectMap<String,Scan> scans=new ObjectMap<String,Scan>(String.class,Scan.class);
	@SuppressWarnings("unchecked")
	final ObjectMap<String,ArrayList<Hit>> partitions=new ObjectMap(String.class,ArrayList.class);
	final ObjectMap<String,long[]> counts=new ObjectMap<String,long[]>(String.class,long[].class);
	final ArrayList<Hit> hits=new ArrayList<Hit>();final ArrayList<Locus> loci=new ArrayList<Locus>();long invalid;
}
