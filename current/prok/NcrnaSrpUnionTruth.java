package prok;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import prok.EukPanelLoci.Hit;
import prok.EukPanelLoci.Locus;
import prok.NcrnaClippingComparison.Contig;
import prok.NcrnaClippingComparison.Feature;
import prok.NcrnaClippingComparison.Genome;
import shared.Shared;
import structures.ByteBuilder;

/** Versioned SRP-only reference overlay using the existing fixed-seed locus rule.
 * Original hits keep their subtype/family bytes; only their grouping is unified.
 * @author Brian, Raiden */
public final class NcrnaSrpUnionTruth {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaSrpUnionTruth.class,false);
		String loci=null,members=null,out=null;long expectedMembers=-1;int expectedOld=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Use flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("members")){members=v;}else if(k.equals("outprefix")){out=v;}
			else if(k.equals("expectedmembers")){expectedMembers=Long.parseLong(v);}else if(k.equals("expectedoldsrp")){expectedOld=Integer.parseInt(v);}
			else{throw new IllegalArgumentException("Unknown argument: "+k);}}
		require(loci!=null && members!=null && out!=null && expectedMembers>=0 && expectedOld>0,"Original loci/members and explicit denominators required");
		final EukPanelLocusDiagnostics original=new EukPanelLocusDiagnostics();original.loadLoci(loci);original.loadMembers(members,expectedMembers);
		final NcrnaSrpUnionTruth x=new NcrnaSrpUnionTruth();final LineParser1 p=new LineParser1('\t');
		for(EukPanelLocusDiagnostics.Locus old:original.loci){if(srp(old.family)){
			x.oldIds.add(old.id+'\t'+old.family);
			for(Hit h:old.members){p.set(h.h.raw);final Hit union=new Hit(h.h.raw,p,"SRP_UNION");
				x.hits.add(union);ArrayList<Hit> group=x.groups.get(union.partition);if(group==null){group=new ArrayList<Hit>();x.groups.put(union.partition,group);}group.add(union);}
		}}
		require(x.oldIds.size()==expectedOld,"Original SRP locus count differs from declared frozen-reference scope");
		x.construct();x.write(out,expectedMembers);Shared.closeStream(pp.outstream);
	}
	void construct(){
		assert(loci.isEmpty() && !hits.isEmpty()):"Rebuild only the declared SRP hits, once";
		final ArrayList<String> keys=new ArrayList<String>();for(String key:groups.keys()){if(key!=null){keys.add(key);}}Collections.sort(keys);
		for(String key:keys){EukPanelLoci.buildPartition(groups.get(key),loci);}
		long n=0;for(Locus locus:loci){EukPanelLoci.findAmbiguity(locus);for(Hit h:locus.members){
			require(h.locus==locus && h.partition.equals(locus.seed.partition) && EukPanelLoci.qualifies(h,locus.seed),"Union membership must obey the existing fixed seed and half-member rule");n++;}}
		require(n==hits.size(),"Every original SRP hit must survive exactly once");
	}
	void write(String prefix,long allMembers){
		assert(!loci.isEmpty()):"Publish only a constructed physical SRP reference";
		RrnaResourceIO.fresh(prefix,".truth.tsv.gz",".members.tsv.gz",".old_ids.tsv",".overlaps.tsv",".summary.tsv");
		final ByteStreamWriter[] writers=new ByteStreamWriter[5];final ByteBuilder b=new ByteBuilder();
		long admitted=0,flagged=0,overlaps=0;
		try{
			final ByteStreamWriter truth=writers[0]=RrnaResourceIO.writer(prefix+".truth.tsv.gz"),members=writers[1]=RrnaResourceIO.writer(prefix+".members.tsv.gz"),old=writers[2]=RrnaResourceIO.writer(prefix+".old_ids.tsv"),pairs=writers[3]=RrnaResourceIO.writer(prefix+".overlaps.tsv"),summary=writers[4]=RrnaResourceIO.writer(prefix+".summary.tsv");
			truth.println(HEADER);members.println(EukPanelTruthDisposition.MEMBER_HEADER+"\tunion_locus_id");old.println("old_locus_id\tfold_family");
			pairs.println("locus_a\tlocus_b\tassembly\tnative_seqid\tstrand\tstart_a\tend_a\tstart_b\tend_b\toverlap_bp\tstatus");
			Collections.sort(oldIds);for(String id:oldIds){old.println(id);}
			for(Hit h:hits){members.print(b.clear().append(h.h.raw).tab().append(id(h.locus)).nl());}
			for(Locus locus:loci){final Hit h=locus.seed,subtype=EukPanelLoci.subtype(locus.members);final boolean ga=subtype!=null,ambiguous=h.spanning || locus.witnessA!=null;
				if(ga){admitted++;}if(ambiguous){flagged++;}
				truth.print(b.clear().append(id(locus)).tab().append(h.h.assembly).tab().append(h.h.nativeSeqid).tab().append(h.h.strand).tab().append(h.h.start).tab().append(h.h.stop)
					.append("\tSRP_UNION\t").append(h.h.bin).tab().append(h.h.compartment).tab().append(h.h.model).tab().append(ga?subtype.h.model:"NA").tab().append(locus.members.size())
					.tab().append(ga?1:0).tab().append(h.h.full?1:0).tab().append(h.h.trunc.equals("no")?0:1).tab().append(h.h.edge).tab().append(geometry(h))
					.tab().append(ambiguous?"SPANNING_AMBIGUOUS":"NONE").tab().append(h.h.gaPass?1:0).tab().append(h.h.gaPass?1:0).tab().append(POLICY).nl());
			}
			// Retain residual overlapping seeds for review; never add an unstated containment merge.
			for(int i=0;i<loci.size();i++){final Hit a=loci.get(i).seed;for(int j=i+1;j<loci.size();j++){final Hit z=loci.get(j).seed;
				if(!a.partition.equals(z.partition)){continue;}final long overlap=EukPanelLoci.overlap(a,z);if(overlap<=0){continue;}overlaps++;
				pairs.print(b.clear().append(id(loci.get(i))).tab().append(id(loci.get(j))).tab().append(a.h.assembly).tab().append(a.h.nativeSeqid).tab().append(a.h.strand)
					.tab().append(a.h.start).tab().append(a.h.stop).tab().append(z.h.start).tab().append(z.h.stop).tab().append(overlap).tab().append("RETAINED_REVIEW_REQUIRED").nl());}}
			summary.println("old_srp_loci\tunion_loci\tsrp_members\tall_original_members\tadmitted\tspanning_flagged\tresidual_overlaps\tpolicy");
			summary.print(b.clear().append(oldIds.size()).tab().append(loci.size()).tab().append(hits.size()).tab().append(allMembers).tab().append(admitted).tab().append(flagged).tab().append(overlaps).tab().append(POLICY).nl());
		}finally{boolean error=false;for(ByteStreamWriter w:writers){if(w!=null){error|=w.poisonAndWait();}}require(!error,"SRP overlay output failed");}
		System.err.println("SRP_UNION_COMPLETE old_loci="+oldIds.size()+" union_loci="+loci.size()+" members="+hits.size()+" flagged="+flagged+" residual_overlaps="+overlaps);
	}
	/** Apply an explicitly bound overlay after frozen-v1 loading, before history/masks. */
	public static void apply(NcrnaClippingComparison grade,String path,String oldPath,int expectedOld,int expectedNew){
		assert(grade.frozenMode && !grade.historyLoaded):"Replace only v1 SRP reference rows before historical masks are loaded";
		final ObjectMap<String,String> old=new ObjectMap<String,String>(String.class,String.class);final LineParser1 p=new LineParser1('\t');
		ByteFile in=RrnaResourceIO.open(oldPath);try{RrnaResourceIO.header(in,"old_locus_id\tfold_family",oldPath);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==2 && srp(p.parseString(1)) && old.put("LOCUS:"+p.parseString(0),p.parseString(1))==null,"Original SRP IDs must be unique and subtype-specific");}}
		finally{require(!in.close(),"Original SRP IDs read failed");}require(old.size()==expectedOld,"Original SRP denominator mismatch");
		int removed=0;for(Genome g:grade.ordered){g.srpReferenceKind=POLICY;for(Contig c:g.ordered){removed+=remove(c.truth,old);removed+=remove(c.masks,old);}}
		require(removed==expectedOld && old.isEmpty(),"Overlay must replace exactly every original SRP locus and no other family");
		int loaded=0;final ObjectSet<String> seen=new ObjectSet<String>(String.class);in=RrnaResourceIO.open(path);
		try{RrnaResourceIO.header(in,HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==21 && p.termEquals("SRP_UNION",6) && p.termEquals(POLICY,20),"Explicit reviewed SRP union policy required");
			final Genome g=grade.genomes.get(p.parseString(1));require(g!=null && p.termEquals(g.panel,7),"Union locus must bind its original assembly/domain");
			final Contig c=g.contigs.get(p.parseString(2));require(c!=null && c.compartment.equals("PROKARYOTIC") && p.termEquals(c.compartment,8),"Union locus contig/compartment changed");
			final String seed=p.parseString(9),subtype=p.parseString(10);final boolean ga=EukPanelLocusDiagnostics.bit(p,12);
			final String family=grade.modelFamilies.get(ga?subtype:seed);require(srp(family) && srp(grade.modelFamilies.get(seed)) && (ga || subtype.equals("NA")),"Union subtype must be an explicitly supported SRP model");
			final String id="SRP_UNION:"+p.parseString(0);require(seen.add(id),"Repeated union locus identity");
			final Feature f=new Feature(c,family,ga?subtype:seed,p.parseInt(4),p.parseInt(5),p.parseString(3),id,p.parseString(17),POLICY);
			require(f.start>0 && f.stop>=f.start && f.stop<=c.length && (f.strand.equals("+") || f.strand.equals("-")),"Union interval must remain inside its physical contig");
			f.geometry=p.parseString(16);require(f.geometry.equals("FULLSPAN_INTACT") || f.geometry.equals("PARTIAL_INTERNAL") || f.geometry.equals("EDGE"),"Unsupported SRP geometry");
			f.truncated=EukPanelLocusDiagnostics.bit(p,14);f.boundaryEligible=EukPanelLocusDiagnostics.bit(p,18);f.bothEligible=EukPanelLocusDiagnostics.bit(p,19);
			require((ga || !f.boundaryEligible) && f.boundaryEligible==f.bothEligible,"SRP endpoint eligibility cannot exceed supported seed evidence");
			if(ga){c.truth.add(f);}else{f.state="BELOW_GA";c.masks.add(f);}loaded++;
		}}finally{require(!in.close(),"Union truth read failed");}
		require(loaded==expectedNew,"Physical SRP locus denominator mismatch");grade.referenceRows+=loaded-removed;grade.expectedFrozenLoci+=loaded-removed;
	}
	static int remove(ArrayList<Feature> features,ObjectMap<String,String> old){
		assert(old!=null):"Removal is constrained by original, explicitly enumerated SRP identities";
		int n=0;for(int i=features.size()-1;i>=0;i--){final Feature f=features.get(i);if(srp(f.family)){require(f.family.equals(old.remove(f.id)),"Unexpected original SRP identity/subtype");features.remove(i);n++;}}return n;
	}
	static String geometry(Hit h){assert(h.h.valid):"Only valid SRP hits enter physical grouping";return !h.h.trunc.equals("no") || !h.h.edge.equals("INTERNAL")?"EDGE":h.h.full?"FULLSPAN_INTACT":"PARTIAL_INTERNAL";}
	static String id(Locus locus){return Integer.toString(locus.id);}
	static boolean srp(String family){return "SRP_small".equals(family) || "SRP_large".equals(family);}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	public static final String POLICY="FRESH_PROK_SRP_UNION_v2";
	public static final String HEADER="locus_id\tassembly\tnative_seqid\tstrand\tstart\tend\tfold_family\tselection_bin\tactual_compartment\tseed_accession\tsubtype_accession\tmembers\tga_pass_any\tfull_model_span\tinfernal_truncated\tcontig_edge_contact\tgeometry_stratum\tflags\tboundary_metric_eligible\tboth_ends_metric_eligible\tpolicy_version";
	final ArrayList<String> oldIds=new ArrayList<String>();final ArrayList<Hit> hits=new ArrayList<Hit>();final ArrayList<Locus> loci=new ArrayList<Locus>();
	@SuppressWarnings("unchecked") final ObjectMap<String,ArrayList<Hit>> groups=new ObjectMap(String.class,ArrayList.class);
}
