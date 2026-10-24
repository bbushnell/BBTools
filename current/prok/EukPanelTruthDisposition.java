package prok;

import java.util.ArrayList;
import java.util.Collections;
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

/** Applies the approved ends-pattern absorptions without reconstructing other loci.
 * All raw hit bytes survive; only locus membership and review dispositions change.
 * The withdrawn 5.8S homolog rule is deliberately absent.
 * @author Raiden */
public final class EukPanelTruthDisposition {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args,EukPanelTruthDisposition.class,false);
		String candidate=null,members=null,cases=null,out=null,exception=null,policy=LEGACY_POLICY;int expectedLoci=-1,expectedMembers=-1,expectedEnds=-1,expectedAmbiguous=-1;
		for(String arg:pp.args){final int eq=arg.indexOf('=');require(eq>0,"Expected flag=value");final String k=arg.substring(0,eq).toLowerCase(Locale.ROOT),v=arg.substring(eq+1);
			if(k.equals("candidate")){candidate=v;}else if(k.equals("members")){members=v;}else if(k.equals("cases")){cases=v;}
			else if(k.equals("outprefix")){out=v;}else if(k.equals("exception")){exception=v;}
			else if(k.equals("expectedloci")){expectedLoci=Integer.parseInt(v);}else if(k.equals("expectedmembers")){expectedMembers=Integer.parseInt(v);}
			else if(k.equals("expectedends")){expectedEnds=Integer.parseInt(v);}else if(k.equals("expectedambiguous")){expectedAmbiguous=Integer.parseInt(v);}
			else if(k.equals("policy")){policy=v;}else{throw new IllegalArgumentException("Unknown flag: "+k);}}
		require(policy.equals(LEGACY_POLICY) || policy.equals(PROK_POLICY),"Unknown reviewed disposition policy");
		final boolean prok=policy.equals(PROK_POLICY);
		require(candidate!=null && members!=null && cases!=null && out!=null && expectedLoci>0 && expectedMembers>0
			&& (prok?exception==null && expectedAmbiguous>=0 && expectedEnds<0:exception!=null && expectedEnds>=0 && expectedAmbiguous<0),"Explicit inputs, policy-specific cases and conservation counts required");
		Shared.setThreads(1);ReadWrite.setZipThreads(1);
		final EukPanelTruthDisposition x=new EukPanelTruthDisposition();x.policy=policy;x.loadCandidate(candidate,expectedLoci);x.loadMembers(members,expectedMembers);
		x.loadCases(cases,exception,prok?expectedAmbiguous:expectedEnds);if(!prok){x.absorb(exception);}x.write(out,expectedLoci,expectedMembers);Shared.closeStream(pp.outstream);
	}
	void loadCandidate(String path,int expected){
		assert(rows.isEmpty()):"One sealed v1.4 candidate supplies the starting locus universe";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,INPUT_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==33 && p.termEquals("FREEZE_CANDIDATE_NOT_FROZEN",32),"Expected unfrozen33-column candidate");
			final Row r=new Row(line,p);
			if(policy.equals(PROK_POLICY)){require(r.compartment.equals("PROKARYOTIC") && (r.bin.equals("bacteria") || r.bin.equals("archaea")),"Fresh prok policy requires original prok compartment/domain");}
			require(byId.put(r.id,r)==null,"Duplicate candidate locus");rows.add(r);
		}}finally{require(!in.close(),"Candidate read failed");}require(rows.size()==expected,"Starting locus denominator changed");
	}
	void loadMembers(String path,int expected){
		assert(entries.isEmpty()):"Every original hit is assigned to exactly one original locus";
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,MEMBER_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);
			require(p.terms()==43 && p.termEquals("REVIEW_ONLY",39),"Expected original40 hit fields and three membership fields");
			final EukPanelLoci.Hit h=new EukPanelLoci.Hit(line,p);final Row r=byId.get(p.parseString(40));
			require(r!=null && h.h.valid && r.assembly.equals(h.h.assembly) && r.seqid.equals(h.h.nativeSeqid) && r.strand.equals(h.h.strand) && r.family.equals(h.h.family) && r.compartment.equals(h.h.compartment),"Member identity differs from owning locus");
			require(p.termEquals(h.h.full?"true":"false",32) && p.termEquals(h.h.gaPass?"true":"false",35),"Member GA/model span differs from raw metrics");
			require(p.termEquals("SEED",41) || p.termEquals("MEMBER",41),"Input may not already contain disposition roles");
			final Entry e=new Entry(h,r,p.termEquals("SEED",41),EukPanelLocusDiagnostics.bit(p,42));
			require(byHit.put(key(h),e)==null,"Repeated original member identity");r.members.add(e);entries.add(e);
			if(e.wasSeed){require(r.seed==null && h.h.run.equals(r.run) && h.table.equals(r.table) && h.h.lineNumber==r.line && h.h.model.equals(r.seedAccession) && h.h.start==r.start && h.h.stop==r.end && h.bits==r.bits,"Original seed identity/geometry changed");r.seed=e;}
		}}finally{require(!in.close(),"Member read failed");}
		require(entries.size()==expected,"Member denominator changed");
		for(Row r:rows){
			require(r.seed!=null && r.members.size()==r.expected,"Each original locus needs its declared members and one seed");
			final ArrayList<EukPanelLoci.Hit> hits=new ArrayList<EukPanelLoci.Hit>();boolean ga=false;
			for(Entry e:r.members){require(EukPanelLoci.qualifies(e.hit,r.seed.hit),"Original membership must satisfy the fixed seed-overlap rule");hits.add(e.hit);ga|=e.hit.h.gaPass;}
			final EukPanelLoci.Hit subtype=EukPanelLoci.subtype(hits);
			require(ga==r.ga && r.subtype.equals(subtype==null?"NA":subtype.h.model),"Original admission/subtype evidence differs from conserved members");
		}
	}
	void loadCases(String path,String exception,int expected){
		assert(parents.isEmpty()):"Absorption parents come only from the reviewed spanning cases";
		final boolean prok=policy.equals(PROK_POLICY);
		final Row except=prok?null:byId.get(exception);require(prok || except!=null && except.ambiguous,"Explicit unchanged exception must be a flagged existing locus");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);
		try{RrnaResourceIO.header(in,EukPanelTruthCandidate.CASE_HEADER,path);for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18,"Spanning-case schema changed");
			final Row r=byId.get(p.parseString(0));require(r!=null && r.ambiguous,"Case must identify an originally flagged locus");
			require(p.termEquals(r.assembly,2) && p.termEquals(r.seqid,3) && p.termEquals(r.strand,4) && p.termEquals(r.family,5) && p.termEquals(r.run,6) && p.parseLong(7)==r.line && p.parseLong(8)==r.start && p.parseLong(9)==r.end,"Spanning case identity changed");
			final Entry a=witness(r,p,10),b=witness(r,p,14);
			require(EukPanelLoci.substantial(a.hit) && EukPanelLoci.substantial(b.hit) && EukPanelLoci.separate(a.hit,b.hit),"Ends witnesses need original GA/full-model and separation evidence");
			if(prok){seen.add(r.id);continue;}// Preserve flagged prok candidates; euk absorptions are not this policy.
			if(r==except){continue;}
			require(r.seed.hit.h.gaPass && r.seed.hit.h.full && (r.seedAccession.equals("RF02543") || r.seedAccession.equals("RF01960")),"Ends-pattern parent must retain its full supported eukaryotic seed");
			if(seen.add(r.id)){parents.add(r);r.ends=true;}
			r.witnesses.add(a);r.witnesses.add(b);
		}}finally{require(!in.close(),"Spanning cases read failed");}
		if(prok){
			require(seen.size()==expected,"Flagged prok locus denominator changed");
			for(Row r:rows){require(r.ambiguous==seen.contains(r.id),"Every flagged prok locus needs its original case evidence");}return;
		}
		require(parents.size()==expected,"Ends-pattern parent count changed");
		for(Row r:rows){if(r.ambiguous){require(r==except || r.ends,"Unclassified spanning locus");}}
	}
	Entry witness(Row r,LineParser1 p,int c){
		assert(c==10 || c==14):"Witness offsets follow the original case schema";
		final Entry e=byHit.get(r.assembly+'\t'+p.parseString(c)+'\t'+p.parseLong(c+1));
		require(e!=null && e.owner.partition.equals(r.partition) && e.owner.family.equals(r.family) && e.hit.h.start==p.parseLong(c+2) && e.hit.h.stop==p.parseLong(c+3) && EukPanelLoci.overlap(e.hit,r.seed.hit)>0,"Witness must be an exact same-fold hit overlapping the parent seed");return e;
	}
	void absorb(String exception){
		assert(!parents.isEmpty()):"This tool applies explicit absorption decisions, not an unsupervised merge heuristic";
		for(Row parent:parents){for(Entry witness:parent.witnesses){
			final Row child=witness.owner;witness.endsWitness=true;
			if(child==parent){continue;}
			require(!child.ends && !child.id.equals(exception),"Cannot absorb another approved euk seed or the unchanged exception");
			require(child.parent==child || child.parent==parent,"Witness locus claimed by two eukaryotic seeds");child.parent=parent;
		}}
		for(Row r:rows){if(r.parent!=r){
			require(r.parent.parent==r.parent && r.family.equals(r.parent.family) && r.partition.equals(r.parent.partition),"Absorption must be acyclic and stay within assembly/contig/strand/fold");
			r.parent.members.addAll(r.members);r.parent.excludeBoth|=r.excludeBoth;absorbed++;
		}}
		long total=0;final ObjectSet<String> assigned=new ObjectSet<String>(String.class);
		for(Row r:rows){if(r.parent==r){for(Entry e:r.members){require(e.owner.parent==r && assigned.add(key(e.hit)),"Final members must occur exactly once under their surviving parent");total++;}}}
		require(total==entries.size(),"Absorption lost or duplicated members");
	}
	void write(String prefix,int oldCount,int memberCount){
		assert(absorbed>=0):"Output conservation is checked against the original immutable universe";
		RrnaResourceIO.fresh(prefix,".truth_candidate.tsv.gz",".locus_members.tsv.gz",".merge_map.tsv",".counts.tsv",".summary.tsv");
		final ByteStreamWriter out=RrnaResourceIO.writer(prefix+".truth_candidate.tsv.gz"),members=RrnaResourceIO.writer(prefix+".locus_members.tsv.gz"),map=RrnaResourceIO.writer(prefix+".merge_map.tsv"),counts=RrnaResourceIO.writer(prefix+".counts.tsv"),summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final LineParser1 p=new LineParser1('\t');final ByteBuilder b=new ByteBuilder();final ObjectMap<String,long[]> bins=new ObjectMap<String,long[]>(String.class,long[].class);long n=0,memberRows=0,witnessRows=0,excluded=0;
		try{
			out.println(OUTPUT_HEADER);members.println(MEMBER_HEADER);map.println("old_locus_id\tnew_locus_id\tassembly\tnative_seqid\tstrand\tfold_family\told_start\told_end\tnew_start\tnew_end\ttransferred_members\treason");
			for(Row r:rows){if(r.parent!=r){final Row z=r.parent;map.print(b.clear().append(r.id).tab().append(z.id).tab().append(r.assembly).tab().append(r.seqid).tab().append(r.strand).tab().append(r.family).tab().append(r.start).tab().append(r.end).tab().append(z.start).tab().append(z.end).tab().append(r.members.size()).tab().append("ENDS_PATTERN_ABSORBED").nl());continue;}
				p.set(r.raw);b.clear();final ArrayList<String> accessions=new ArrayList<String>();final ObjectSet<String> seen=new ObjectSet<String>(String.class);final ArrayList<EukPanelLoci.Hit> hits=new ArrayList<EukPanelLoci.Hit>();int spanning=0;
				for(Entry e:r.members){hits.add(e.hit);if(seen.add(e.hit.h.model)){accessions.add(e.hit.h.model);}if(e.spanning){spanning++;}}
				Collections.sort(accessions);final EukPanelLoci.Hit subtype=EukPanelLoci.subtype(hits);
				for(int i=0;i<33;i++){if(i>0){b.tab();}
					if(i==18){b.append(r.members.size());}else if(i==19){for(int j=0;j<accessions.size();j++){if(j>0){b.append(';');}b.append(accessions.get(j));}}
					else if(i==21){b.append(spanning);}else if(i==25){b.append(subtype==null?"NA":subtype.h.model);}
					else if(i==28){String flags=r.flags;if(r.ends){flags=flags.replace("SPANNING_AMBIGUOUS","SPANNING_ENDS_PATTERN");}if(r.excludeBoth && !EukPanel58sHomologAudit.member(flags,"AMBIGUOUS_CROSSFOLD")){flags+=";AMBIGUOUS_CROSSFOLD";}b.append(flags);}
					else{p.setBounds(i);b.append(p.line(),p.a(),p.b()-p.a());}
				}
				final boolean both=r.boundaryEligible && !r.excludeBoth;out.print(b.tab().append(both?1:0).tab().append(policy).nl());n++;if(!both){excluded++;}
				EukPanelCmProposal.increment(bins,r.bin+'\t'+r.family+'\t'+r.compartment+'\t'+r.admission+'\t'+r.geometry,1);
			}
			for(Entry e:entries){final Row r=e.owner.parent;p.set(e.hit.h.raw);p.setBounds(39);
				final String role=e==r.seed?"SEED":e.endsWitness?"ENDS_WITNESS":"MEMBER";
				members.print(b.clear().append(e.hit.h.raw,0,p.b()).tab().append(r.id).tab().append(role).tab().append(e.spanning?1:0).nl());memberRows++;if(role.equals("ENDS_WITNESS")){witnessRows++;}
			}
			counts.println("selection_bin\tfold_family\tactual_compartment\tadmission_class\tgeometry_stratum\tloci");long countSum=0;
			for(String k:EukPanelCmProposal.sorted(bins)){final long value=bins.get(k)[0];counts.print(b.clear().append(k).tab().append(value).nl());countSum+=value;}
			require(n==oldCount-absorbed && n==countSum && memberRows==memberCount,"Final locus/member denominator algebra failed");
			summary.println("original_loci\tsurviving_loci\tabsorbed_loci\tmembers\tends_pattern_parents\tends_witness_rows\tboth_end_ineligible_loci\tstatus");
			summary.print(b.clear().append(oldCount).tab().append(n).tab().append(absorbed).tab().append(memberRows).tab().append(parents.size()).tab().append(witnessRows).tab().append(excluded).tab().append("FREEZE_CANDIDATE_NOT_FROZEN").nl());
		}finally{RrnaResourceIO.finish(out,members,map,counts,summary);}
		System.err.println("TRUTH_DISPOSITION_COMPLETE original_loci="+oldCount+" loci="+n+" absorbed="+absorbed+" members="+memberRows+" ends_parents="+parents.size()+" frozen=false");
	}
	static String key(EukPanelLoci.Hit h){return h.h.assembly+'\t'+h.h.run+'\t'+h.h.lineNumber;}
	static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	static final class Row{
		Row(byte[] line,LineParser1 p){raw=line;id=p.parseString(0);assembly=p.parseString(1);seqid=p.parseString(2);strand=p.parseString(3);start=p.parseLong(4);end=p.parseLong(5);family=p.parseString(6);bin=p.parseString(7);compartment=p.parseString(8);seedAccession=p.parseString(9);run=p.parseString(10);table=p.parseString(11);this.line=p.parseLong(12);bits=p.parseDouble(13);ga=EukPanelLocusDiagnostics.bit(p,14);expected=p.parseInt(18);subtype=p.parseString(25);admission=p.parseString(26);geometry=p.parseString(27);flags=p.parseString(28);boundaryEligible=EukPanelLocusDiagnostics.bit(p,31);ambiguous=EukPanel58sHomologAudit.member(flags,"SPANNING_AMBIGUOUS");excludeBoth=EukPanel58sHomologAudit.member(flags,"AMBIGUOUS_CROSSFOLD");partition=assembly+'\t'+seqid+'\t'+strand;parent=this;
			require(start>0 && end>=start && expected>0 && (strand.equals("+") || strand.equals("-")),"Invalid original locus geometry");require(admission.equals(ga?"ADMITTED":"BELOW_GA"),"Original candidate admission must equal any-GA policy");}
		final byte[] raw;final String id,assembly,seqid,strand,family,bin,compartment,seedAccession,run,table,subtype,admission,geometry,flags,partition;final long start,end,line;final double bits;final int expected;final boolean ga,ambiguous,boundaryEligible;
		Row parent;Entry seed;boolean ends,excludeBoth;final ArrayList<Entry> members=new ArrayList<Entry>(),witnesses=new ArrayList<Entry>();
	}
	static final class Entry{Entry(EukPanelLoci.Hit h,Row r,boolean seed,boolean span){hit=h;owner=r;wasSeed=seed;spanning=span;}final EukPanelLoci.Hit hit;final Row owner;final boolean wasSeed,spanning;boolean endsWitness;}
	public static final String INPUT_HEADER=EukPanelTruthCandidate.CANDIDATE_HEADER;
	public static final String OUTPUT_HEADER=INPUT_HEADER+"\tboth_ends_metric_eligible\tpolicy_version";
	public static final String MEMBER_HEADER=EukPanelLoci.HEADER+"\tlocus_id\trole\tspanning_hit";
	public static final String CASE_HEADER=EukPanelTruthCandidate.CASE_HEADER;
	public static final String LEGACY_POLICY="1.5_WITH_HOMOLOG_RULE_WITHDRAWN",PROK_POLICY="FRESH_PROK_ADMISSION_v1";
	String policy=LEGACY_POLICY;
	final ArrayList<Row> rows=new ArrayList<Row>(),parents=new ArrayList<Row>();final ArrayList<Entry> entries=new ArrayList<Entry>();
	final ObjectMap<String,Row> byId=new ObjectMap<String,Row>(String.class,Row.class);final ObjectMap<String,Entry> byHit=new ObjectMap<String,Entry>(String.class,Entry.class);int absorbed;
}
