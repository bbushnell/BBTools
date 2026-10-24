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

/** Q4 review cases from descriptive loci and their complete original members.
 * No admission or scoring. Crossfold cases require two GA-supported SSU/LSU
 * seed intervals with positive same-strand overlap. Unexpected compartment is
 * the draft's declared NUCLEAR RF00010/RF00023 case definition, not a veto.
 * @author Ganyu */
public final class EukPanelLocusDiagnostics {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelLocusDiagnostics.class,false);
		String loci=null,members=null,out=null;long expected=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");
			final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("loci")){loci=v;}else if(k.equals("members")){members=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedmembers")){expected=Long.parseLong(v);}
			else{throw new IllegalArgumentException("Unknown flag: "+k);}
		}
		require(loci!=null && members!=null && out!=null && expected>=0,"Explicit loci, members, output and conserved denominator required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelLocusDiagnostics x=new EukPanelLocusDiagnostics();x.loadLoci(loci);x.loadMembers(members,expected);x.write(out);Shared.closeStream(pp.outstream);
	}
	void loadLoci(String path){
		assert(byId.isEmpty()):"One descriptive locus generation supplies all diagnostic keys";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,LOCUS_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==26 && p.termEquals("DESCRIPTIVE_ONLY",23),"Expected descriptive 26-column locus table with distinct geometry and subtype");
			final Locus x=new Locus(p);require(byId.put(x.id,x)==null,"Duplicate locus ID: "+x.id);loci.add(x);
		}}finally{require(!in.close(),"Locus read failed");}
	}
	void loadMembers(String path,long expected){
		assert(memberCount==0):"Each original member must be checked once against its locus";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> identities=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,EukPanelLoci.HEADER+"\tlocus_id\trole\tspanning_hit",path);
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);require(p.terms()==43 && p.termEquals("REVIEW_ONLY",39),"Expected original 40-column hit plus three membership columns");
				final EukPanelLoci.Hit h=new EukPanelLoci.Hit(line,p);final Locus x=byId.get(p.parseString(40));
				require(x!=null && h.h.valid && x.assembly.equals(h.h.assembly) && x.seqid.equals(h.h.nativeSeqid) && x.strand.equals(h.h.strand)
					&& x.family.equals(h.h.family) && x.compartment.equals(h.h.compartment),"Member identity/compartment disagrees with locus");
				require(identities.add(h.h.assembly+'\t'+h.h.run+'\t'+h.table+'\t'+h.h.lineNumber),"Duplicate original member identity");
				require(p.termEquals(h.h.full?"true":"false",32) && p.termEquals(h.h.gaPass?"true":"false",35),"Derived member flags disagree with original metrics");
				require(p.termEquals("SEED",41) || p.termEquals("MEMBER",41),"Unknown member role");
				if(p.termEquals("SEED",41)){
					require(x.seed==null && x.run.equals(h.h.run) && x.table.equals(h.table) && x.line==h.h.lineNumber && x.start==h.h.start && x.end==h.h.stop
						&& x.accession.equals(h.h.model) && x.seedBits==h.bits && x.seedGa==h.h.gaPass,"Seed identity or metrics disagree with original member");x.seed=h;
				}
				x.members.add(h);x.observedGa|=h.h.gaPass;memberCount++;
			}
		}finally{require(!in.close(),"Member read failed");}
		require(memberCount==expected,"Member denominator changed: "+memberCount+" != "+expected);
		for(Locus x:loci){
			require(x.seed!=null && x.members.size()==x.expected && x.observedGa==x.ga,"Locus member conservation or any-GA flag disagrees: "+x.id);
			final EukPanelLoci.Hit subtype=EukPanelLoci.subtype(x.members);
			require(x.subtype.equals(subtype==null?"NA":subtype.h.model),"Subtype accession must equal highest-bit GA member independently of geometry: "+x.id);
			Collections.sort(x.members,EukPanelLoci.RANK);
		}
	}
	void write(String prefix){
		assert(memberCount>=0):"Diagnostics preserve the validated member denominator; they never admit truth";
		RrnaResourceIO.fresh(prefix,".diagnostic_cases.tsv",".diagnostic_counts.tsv");
		final ByteStreamWriter cases=RrnaResourceIO.writer(prefix+".diagnostic_cases.tsv"),counts=RrnaResourceIO.writer(prefix+".diagnostic_counts.tsv");
		final ByteBuilder b=new ByteBuilder();long seedBelow=0,unexpected=0,unexpectedRows=0,pairs=0;
		final ObjectSet<String> crossLoci=new ObjectSet<String>(String.class);
		Collections.sort(loci,ORDER);
		try{
			cases.println("reason\tlocus_a\tlocus_b\tassembly\tnative_seqid\tstrand\tcompartment\toverlap_bp\tfamily_a\tmodel_a\tstart_a\tend_a\trun_a\ttable_a\tline_a\tbits_a\tga_a\tmodel_from_a\tmodel_to_a\tclen_a\ttrunc_a\tfamily_b\tmodel_b\tstart_b\tend_b\trun_b\ttable_b\tline_b\tbits_b\tga_b\tmodel_from_b\tmodel_to_b\tclen_b\ttrunc_b\tstatus");
			for(int i=0;i<loci.size();i++){
				final Locus x=loci.get(i);if(!x.ga){continue;}
				if(!x.seedGa){
					EukPanelLoci.Hit supporting=null;for(EukPanelLoci.Hit h:x.members){if(h.h.gaPass){supporting=h;break;}}
					require(supporting!=null,"Any-GA locus needs an actual supporting member");
					writeCase(cases,b,"SEED_BELOW_GA",x,x,x.seed,supporting);seedBelow++;
				}
				boolean found=false;
				if(x.compartment.equals("NUCLEAR")){for(EukPanelLoci.Hit h:x.members){
					if(h.h.gaPass && (h.h.model.equals("RF00010") || h.h.model.equals("RF00023"))){
						writeCase(cases,b,"COMPARTMENT_UNEXPECTED",x,null,h,null);found=true;unexpectedRows++;
					}
				}}if(found){unexpected++;}
				if(!ribosomal(x)){continue;}
				for(int j=i+1;j<loci.size();j++){
					final Locus y=loci.get(j);if(!x.partition.equals(y.partition) || y.start>x.end){break;}
					if(y.ga && ribosomal(y) && !x.family.equals(y.family)){
						writeCase(cases,b,"AMBIGUOUS_CROSSFOLD",x,y,x.seed,y.seed);pairs++;crossLoci.add(x.id);crossLoci.add(y.id);
					}
				}
			}
			counts.println("loci\tmembers\tseed_below_ga_loci\tunexpected_compartment_loci\tunexpected_supporting_hits\tcrossfold_pairs\tcrossfold_loci");
			counts.print(b.clear().append(loci.size()).tab().append(memberCount).tab().append(seedBelow).tab().append(unexpected).tab().append(unexpectedRows).tab().append(pairs).tab().append(crossLoci.size()).nl());
		}finally{RrnaResourceIO.finish(cases,counts);}
		System.err.println("LOCUS_DIAGNOSTICS_COMPLETE loci="+loci.size()+" members="+memberCount+" seed_below_ga="+seedBelow+" unexpected_loci="+unexpected+" crossfold_pairs="+pairs+" no_truth_admission=true");
	}
	static void writeCase(ByteStreamWriter out,ByteBuilder b,String reason,Locus x,Locus y,EukPanelLoci.Hit a,EukPanelLoci.Hit z){
		assert(x!=null && a!=null):"Every case has an original hit witness and a validated locus";
		b.clear().append(reason).tab().append(x.id).tab().append(y==null?"NA":y.id).tab().append(x.assembly).tab().append(x.seqid).tab().append(x.strand).tab().append(x.compartment)
			.tab().append(z==null?0:EukPanelLoci.overlap(a,z));appendHit(b,a);appendHit(b,z);out.print(b.tab().append("REVIEW_ONLY").nl());
	}
	static void appendHit(ByteBuilder b,EukPanelLoci.Hit h){
		assert(b!=null):"Case fields append to the same reusable row buffer";
		if(h==null){for(int i=0;i<13;i++){b.tab().append("NA");}return;}
		final LineParser1 p=new LineParser1('\t');p.set(h.h.raw);
		b.tab().append(h.h.family).tab().append(h.h.model).tab().append(h.h.start).tab().append(h.h.stop).tab().append(h.h.run).tab().append(h.table).tab().append(h.h.lineNumber)
			.tab().append(h.bits,6).tab().append(p.parseDouble(26),6).tab().append(h.h.mf).tab().append(h.h.mt).tab().append(h.h.clen).tab().append(h.h.trunc);
	}
	static boolean ribosomal(Locus x){return x.family.equals("SSU") || x.family.equals("LSU");}
	static boolean bit(LineParser1 p,int i){require(p.termEquals("0",i) || p.termEquals("1",i),"Expected 0/1 locus flag");return p.termEquals("1",i);}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Locus{
		Locus(LineParser1 p){
			id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);start=p.parseLong(4);end=p.parseLong(5);
			family=p.parseString(6);compartment=p.parseString(8);accession=p.parseString(9);run=p.parseString(10);table=p.parseString(11);line=p.parseLong(12);
			seedBits=p.parseDouble(13);ga=bit(p,14);seedGa=bit(p,24);subtype=p.parseString(25);expected=p.parseInt(18);partition=assembly+'\t'+seqid+'\t'+strand;
			require(!id.isEmpty() && start>0 && end>=start && expected>0 && (strand.equals("+") || strand.equals("-")) && Double.isFinite(seedBits),"Invalid descriptive locus");
		}
		final String id,assembly,seqid,strand,family,compartment,accession,subtype,run,table,partition;final long start,end,line;final double seedBits;final boolean ga,seedGa;final int expected;
		boolean observedGa;EukPanelLoci.Hit seed;final ArrayList<EukPanelLoci.Hit> members=new ArrayList<EukPanelLoci.Hit>();
	}
	static final Comparator<Locus> ORDER=(a,b)->{int c=a.partition.compareTo(b.partition);if(c==0){c=Long.compare(a.start,b.start);}if(c==0){c=Long.compare(a.end,b.end);}if(c==0){c=a.id.compareTo(b.id);}return c;};
	public static final String LOCUS_HEADER=EukPanelLoci.LOCUS_HEADER;
	final ObjectMap<String,Locus> byId=new ObjectMap<String,Locus>(String.class,Locus.class);final ArrayList<Locus> loci=new ArrayList<Locus>();long memberCount;
}
