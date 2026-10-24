package prok;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Reconcile saved inclusive P1 SSU/LSU truth rows with exact historical complete
 * membership and original GFF metadata. No calling, grading or truth filtering.
 * @author Ganyu */
public final class NcrnaP1MissInventory {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaP1MissInventory.class,false);
		String refs=null,audit=null,ssu=null,lsu=null,out=null;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(java.util.Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("refs")){refs=v;}else if(k.equals("audit")){audit=v;}else if(k.equals("complete16")){ssu=v;}else if(k.equals("complete23")){lsu=v;}else if(k.equals("outprefix")){out=v;}else{throw new IllegalArgumentException(k);}
		}
		require(refs!=null && audit!=null && ssu!=null && lsu!=null && out!=null,"Explicit frozen references, audit and complete views required");
		Shared.setThreads(1);fileIO.ReadWrite.setZipThreads(1);
		final NcrnaP1MissInventory x=new NcrnaP1MissInventory();x.loadComplete(ssu,"16S",290);x.loadComplete(lsu,"23S",288);x.loadSources(refs);x.run(audit,out,5,383,370);Shared.closeStream(pp.outstream);
	}
	public void loadComplete(String path,String family,int expected){
		assert(expected>0):"The exact sealed V10 reference count is an input contract, not a newly chosen filter";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');int n=0;
		try{RrnaResourceIO.header(in,"seqid\tfamily\tstart\tstop\tstrand",path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==5 && p.termEquals(family,1),"Complete-reference schema/family changed");
			require(complete.add(key(p.parseString(0),NcrnaClippingComparison.fold(family),p.parseInt(2),p.parseInt(3),p.parseString(4))),"Duplicate complete-reference key");n++;
		}}finally{require(!in.close(),"Complete-reference read failed");}require(n==expected,"Complete-reference denominator changed");
	}
	void loadSources(String refs){
		assert(sources.isEmpty()):"Source IDs use immutable GFF path and original physical line number";
		final ByteFile in=RrnaResourceIO.open(refs);final LineParser1 p=new LineParser1('\t');int n=0;
		try{RrnaResourceIO.header(in,"index\tgenome\tfasta\tfasta_sha80\tgff\tgff_sha80",refs);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==6,"Reference manifest schema changed");loadGff(p.parseString(4),p.parseString(1));n++;}}
		finally{require(!in.close(),"Reference manifest read failed");}require(n==203,"Expected original203 source manifests");
	}
	public void loadGff(String path,String assembly){
		assert(path!=null && assembly!=null):"Every source annotation remains attributable to its exact assembly";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');long number=0;
		try{for(byte[] line;(line=in.nextLine())!=null;){number++;if(line.length==0 || line[0]=='#'){continue;}p.set(line);require(p.terms()==9,"Source GFF requires9 fields");
			if(!NcrnaClippingComparison.rna(p)){continue;}final Source s=new Source(assembly,p.parseString(0),p.parseInt(3),p.parseInt(4),p.parseString(6),p.parseString(8));
			require(sources.put(path+":"+number,s)==null,"Duplicate original source identity");
		}}finally{require(!in.close(),"Source GFF read failed");}
	}
	public void run(String audit,String prefix,int expectedArms,int expectedSsu,int expectedLsu){
		assert(expectedArms>0 && !sources.isEmpty()):"Inventory consumes already graded source rows and never invents missing annotations";
		RrnaResourceIO.fresh(prefix,".misses.tsv",".summary.tsv");
		final ByteStreamWriter misses=RrnaResourceIO.writer(prefix+".misses.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteFile in=RrnaResourceIO.open(audit);final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();
		final ObjectSet<String> seen=new ObjectSet<String>(String.class),arms=new ObjectSet<String>(String.class),observedComplete=new ObjectSet<String>(String.class);
		final ObjectMap<String,long[]> counts=new ObjectMap<String,long[]>(String.class,long[].class),totals=new ObjectMap<String,long[]>(String.class,long[].class);
		try{RrnaResourceIO.header(in,AUDIT_HEADER,audit);
			misses.println("arm\tassembly\tcontig\tfold\tsource_id\tstart\tstop\tstrand\tsource_length\tgeometry\texact_V10_complete_member\tsource_partial_attribute\toriginal_attributes");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18,"Saved audit schema changed");
				if(!p.termEquals("P1",1) || !p.termEquals("TRUTH",6) || !(p.termEquals("SSU",5) || p.termEquals("LSU",5))){continue;}
				final String arm=p.parseString(0),family=p.parseString(5),id=p.parseString(7),geometry=p.parseString(16);final Source s=sources.get(id);
				require(s!=null && p.termEquals(s.assembly,2) && p.termEquals(s.seqid,3) && p.parseInt(9)==s.start && p.parseInt(10)==s.stop && p.termEquals(s.strand,11),"Audit/source coordinates or identity disagree: "+id);
				require(p.termEquals("SOURCE_ANNOTATION_RELATIVE",17) && (p.termEquals("TP",12) || p.termEquals("FN",12)),"Expected original inclusive source truth outcomes");
				require(seen.add(arm+'\t'+id),"Repeated source truth in arm");arms.add(arm);
				final String key=key(s.seqid,family,s.start,s.stop,s.strand);final boolean member=complete.contains(key),miss=p.termEquals("FN",12);if(member){observedComplete.add(arm+'\t'+key);}
				String partial=EukPanelCmProposal.gffAttribute(s.attrs,"partial");if(partial==null){partial="UNSPECIFIED";}
				final String group=arm+'\t'+family+'\t'+geometry+'\t'+(member?"YES":"NO")+'\t'+partial;
				long[] v=counts.get(group);if(v==null){v=new long[2];counts.put(group,v);}v[miss?1:0]++;EukPanelCmProposal.increment(totals,arm+'\t'+family,1);
				if(miss){misses.print(b.clear().append(arm).tab().append(s.assembly).tab().append(s.seqid).tab().append(family).tab().append(id).tab().append(s.start).tab().append(s.stop).tab().append(s.strand)
					.tab().append((long)s.stop-s.start+1).tab().append(geometry).tab().append(member?"YES":"NO").tab().append(partial).tab().append(s.attrs).nl());}
			}
			require(arms.size()==expectedArms && observedComplete.size()==(long)expectedArms*complete.size(),"Expected arms and all exact complete keys must appear in each arm source truth");
			for(String arm:arms.keys()){if(arm!=null){require(EukPanelCmProposal.get(totals,arm+"\tSSU")==expectedSsu && EukPanelCmProposal.get(totals,arm+"\tLSU")==expectedLsu,"Inclusive source denominator changed: "+arm);}}
			summary.println("arm\tfold\tgeometry\texact_V10_complete_member\tsource_partial_attribute\tTP\tFN");
			for(String key:EukPanelCmProposal.sorted(counts)){final long[] v=counts.get(key);summary.print(b.clear().append(key).tab().append(v[0]).tab().append(v[1]).nl());}
		}finally{final boolean error=in.close();RrnaResourceIO.finish(misses,summary);require(!error,"Audit read failed");}
		System.err.println("P1_MISS_INVENTORY_COMPLETE arms="+arms.size()+" ssu_per_arm="+expectedSsu+" lsu_per_arm="+expectedLsu+" no_new_grading=true");
	}
	static String key(String seq,String family,int start,int stop,String strand){return seq+'\t'+family+'\t'+start+'\t'+stop+'\t'+strand;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Source{
		Source(String a,String q,int l,int r,String s,String at){assembly=a;seqid=q;start=l;stop=r;strand=s;attrs=at;}
		final String assembly,seqid,strand,attrs;final int start,stop;
	}
	public static final String AUDIT_HEADER="arm\tpanel\tassembly\tcontig\tcompartment\tview\trole\tid\tlabel\tstart\tstop\tstrand\toutcome\tpartner\td5\td3\tgeometry\treference_kind";
	final ObjectMap<String,Source> sources=new ObjectMap<String,Source>(String.class,Source.class);
	final ObjectSet<String> complete=new ObjectSet<String>(String.class);
}
