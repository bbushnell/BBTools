package prok;

import java.util.ArrayList;
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

/** Measures the v1.5 homolog predicate before any truth mutation.
 * Merging only LSU witness loci into RF02543 seeds cannot increase this count:
 * those merged loci have an euk-LSU member and are explicitly ineligible.
 * SSU merges do not change LSU or 5.8S eligibility. Therefore a below-range count
 * is already a failing upper bound; an in-range count is NOT freeze approval.
 * @author Raiden */
public final class EukPanel58sHomologAudit {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanel58sHomologAudit.class,false);
		String in=null,out=null;int expected=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("in")){in=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("expectedloci")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(in!=null && out!=null && expected>=0,"Explicit candidate and denominator required");Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanel58sHomologAudit x=new EukPanel58sHomologAudit();x.load(in,expected);x.write(out);Shared.closeStream(pp.outstream);
	}
	void load(String path,int expected){
		assert(rows.isEmpty()):"Only one immutable candidate supplies the audit universe";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,EukPanelTruthCandidate.CANDIDATE_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){
			p.set(line);require(p.terms()==33 && p.termEquals("FREEZE_CANDIDATE_NOT_FROZEN",32),"Expected unfrozen v1.4 candidate schema");
			final Ref r=new Ref(p);require(seen.add(r.id),"Duplicate locus id");rows.add(r);
			long[] c=counts.get(r.bin);if(c==null){c=new long[5];counts.put(r.bin,c);}
			if(r.family.equals("LSU")){c[2]++;if(eligible(r)){c[3]++;}}
			if(r.family.equals("5.8S") && r.admitted){c[0]++;if(r.crossFamily){c[1]++;}
				ArrayList<Ref> list=r58.get(r.partition);if(list==null){list=new ArrayList<Ref>();r58.put(r.partition,list);}list.add(r);}
		}}finally{require(!in.close(),"Candidate read failed");}require(rows.size()==expected,"Candidate locus denominator changed");
	}
	void write(String prefix){
		assert(!rows.isEmpty()):"The diagnostic is an upper-bound calculation over existing loci, not an admission operation";
		RrnaResourceIO.fresh(prefix,".counts.tsv",".matches.tsv");final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".counts.tsv"),matches=RrnaResourceIO.writer(prefix+".matches.tsv");final ByteBuilder b=new ByteBuilder();long matched=0;
		try{matches.println("selection_bin\tcompartment\tassembly\tnative_seqid\tstrand\tlsu_locus\tlsu_subtype\tlsu_start\tlsu_end\tr58_locus\tr58_start\tr58_end\toverlap_bp\tstatus");
			for(Ref r:rows){if(!eligible(r)){continue;}final ArrayList<Ref> list=r58.get(r.partition);if(list==null){continue;}boolean found=false;
				for(Ref z:list){if(qualifies(r,z)){
					matches.print(b.clear().append(r.bin).tab().append(r.compartment).tab().append(r.assembly).tab().append(r.seqid).tab().append(r.strand).tab().append(r.id)
						.tab().append(r.subtype).tab().append(r.start).tab().append(r.end).tab().append(z.id).tab().append(z.start).tab().append(z.end).tab().append(overlap(r,z)).tab().append("AUDIT_ONLY_NO_TRUTH_CHANGE").nl());found=true;}}
				if(found){counts.get(r.bin)[4]++;matched++;}
			}
			out.println("selection_bin\tadmitted_r58_loci\tr58_cross_family_loci\tlsu_loci\teligible_no_euk_member_lsu\tmatched_lsu_upper_bound\tminimum_85percent\tmaximum_115percent\tassessment");
			for(String key:EukPanelCmProposal.sorted(counts)){final long[] c=counts.get(key);final long lower=(85*c[1]+99)/100,upper=115*c[1]/100;
				b.clear().append(key);for(long n:c){b.tab().append(n);}out.print(b.tab().append(lower).tab().append(upper).tab().append(c[4]<lower?"BELOW_RANGE_BEFORE_MERGE":c[4]>upper?"ABOVE_RANGE_REQUIRES_MERGED_COUNT":"IN_RANGE_STILL_REQUIRES_MERGED_COUNT").nl());}
		}finally{RrnaResourceIO.finish(out,matches);}
		System.err.println("HOMOLOG_AUDIT_COMPLETE loci="+rows.size()+" matched_lsu_upper_bound="+matched+" no_truth_change=true");
	}
	static boolean eligible(Ref r){return r.family.equals("LSU") && !r.eukMember && (r.subtype.equals("RF02541") || r.subtype.equals("RF02540"));}
	static boolean qualifies(Ref a,Ref b){
		assert(a.start>0 && a.end>=a.start && b.start>0 && b.end>=b.start):"Reciprocal overlap uses positive inclusive source coordinates";
		return b.family.equals("5.8S") && b.admitted && a.partition.equals(b.partition) && overlap(a,b)>=(a.end-a.start+2)/2 && overlap(a,b)>=(b.end-b.start+2)/2;
	}
	static long overlap(Ref a,Ref b){return Math.max(0,Math.min(a.end,b.end)-Math.max(a.start,b.start)+1);}
	static boolean member(String list,String accession){for(String s:list.split(";")){if(s.equals(accession)){return true;}}return false;}
	public static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Ref{
		Ref(LineParser1 p){id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);start=p.parseLong(4);end=p.parseLong(5);family=p.parseString(6);bin=p.parseString(7);compartment=p.parseString(8);subtype=p.parseString(25);
			partition=assembly+'\t'+seqid+'\t'+strand;admitted=p.termEquals("ADMITTED",26);eukMember=member(p.parseString(19),"RF02543") || member(p.parseString(19),"RF02546");crossFamily=member(p.parseString(28),"CROSS_FAMILY_OVERLAP");
			require(start>0 && end>=start && (strand.equals("+") || strand.equals("-")),"Invalid candidate interval");require(admitted || p.termEquals("BELOW_GA",26) || p.termEquals("OUT_OF_MAP",26),"Unexpected admission class");
		}
		final String id,assembly,seqid,strand,family,bin,compartment,subtype,partition;final long start,end;final boolean admitted,eukMember,crossFamily;
	}
	final ArrayList<Ref> rows=new ArrayList<Ref>();final ObjectMap<String,long[]> counts=new ObjectMap<String,long[]>(String.class,long[].class);
	@SuppressWarnings("unchecked")
	final ObjectMap<String,ArrayList<Ref>> r58=new ObjectMap(String.class,ArrayList.class);
}
