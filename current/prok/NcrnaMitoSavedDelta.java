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
import parse.Parse;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Joins saved native audit identities, GFF deltas and earlier origin evidence.
 * No alignment, rematching, truth edit or new origin inference.
 * @author Raiden */
public final class NcrnaMitoSavedDelta {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaMitoSavedDelta.class, false);
		final ObjectMap<String,String> f=new ObjectMap<String,String>(String.class, String.class);
		try{
			for(String arg:pp.args){final String[] p=Parse.splitOnFirst(arg, '=');require(p.length==2, "Expected flag=value");
				final String k=p[0].toLowerCase(Locale.ROOT);require(KNOWN.contains("|"+k+"|") && f.put(k, p[1])==null, "Unknown/repeated argument: "+k);}
			for(String k:REQUIRED){require(f.get(k)!=null, "Missing saved input: "+k);}
			Shared.setThreads(1);ReadWrite.setZipThreads(1);
			final NcrnaMitoSavedDelta x=new NcrnaMitoSavedDelta();
			x.audit(f.get("audit"));x.summary(f.get("summary"));x.old(f.get("old"));
			x.delta(f.get("p2delta"), "P2");x.delta(f.get("p3delta"), "P3");x.write(f.get("outprefix"));
		}finally{Shared.closeStream(pp.outstream);}
	}

	public void audit(String path){
		require(refs.isEmpty(), "One immutable audit per comparison");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in, NcrnaEndpointAudit.HEADER, path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==18, "Native audit must retain18 columns");
				if(!(p.termEquals("P2",1) || p.termEquals("P3",1)) || !(p.termEquals("SSU",5) || p.termEquals("tRNA",5))){continue;}
				final int arm=p.termEquals("baseline",0)?0:p.termEquals("C",0)?1:-1;require(arm>=0, "Only the selected baseline/C comparison is supported");
				final Row r=new Row(p);Counts c=counts.get(r.group());if(c==null){c=new Counts();counts.put(r.group(),c);}
				if(p.termEquals("TRUTH",6)){
					require(r.outcome.equals("TP") || r.outcome.equals("FN"), "Reference outcome must be TP/FN");
					final String key=r.group()+'\t'+r.assembly+'\t'+r.contig+'\t'+r.id;
					Row[] pair=refs.get(key);if(pair==null){pair=new Row[2];refs.put(key,pair);}
					require(pair[arm]==null, "Duplicate reference identity within arm");pair[arm]=r;c.n[arm][0]++;
					if(r.outcome.equals("TP")){c.n[arm][1]++;require(!r.partner.equals("NA"),"TP requires its original matched call ID");}
				}else{
					require(p.termEquals("CALL",6), "Unknown audit role");require(calls.put(r.id,r)==null, "Original call IDs must be unique in SSU/tRNA views");
					if(arm==0){baseline.add(r.coordinate());if(r.fp()){oldFp.add(r.coordinate());}}
					else{candidate.add(r.coordinate());if(r.fp() && r.family.equals("SSU")){unmatched.add(r);}}
					if(r.fp()){c.n[arm][2]++;}
				}
			}
		}finally{require(!in.close(), "Read the complete saved audit");}
		for(Row[] pair:refs.values()){if(pair!=null){final Row a=pair[0],z=pair[1];
			if(a==null || z==null){throw new IllegalArgumentException("Missing arm for reference");}
			require(a.coordinate().equals(z.coordinate()) && a.group().equals(z.group()) && a.id.equals(z.id), "Reference identity/geometry must be immutable");
			final Counts c=counts.get(a.group());if(a.outcome.equals("TP") && z.outcome.equals("FN")){c.lost++;losses.add(a);}
			if(a.outcome.equals("FN") && z.outcome.equals("TP")){c.gained++;}
		}}
		require(!refs.isEmpty(), "Empty selected audit population");
	}
	public void summary(String path){
		require(!counts.isEmpty(), "Audit precedes independent summary reconciliation");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');boolean header=false;
		try{for(byte[] line;(line=in.nextLine())!=null;){if(line.length>0 && line[0]=='#'){continue;}p.set(line);
			if(!header){require(p.terms()==19 && p.termEquals("TP",4) && p.termEquals("apparentFP",5), "Expected native19-column summary");header=true;continue;}
			require(p.terms()==19,"Truncated summary");final String key=p.parseString(1)+'\t'+p.parseString(2)+'\t'+p.parseString(3);final Counts c=counts.get(key);if(c==null){continue;}
			final int arm=p.termEquals("baseline",0)?0:p.termEquals("C",0)?1:-1;require(arm>=0 && !c.checked[arm],"Duplicate/unknown summary arm");
			require(c.n[arm][1]==p.parseLong(4) && c.n[arm][2]==p.parseLong(5) && c.n[arm][0]-c.n[arm][1]==p.parseLong(6),"Audit TP/FP/FN differs from summary: "+key);c.checked[arm]=true;
		}}finally{require(!in.close(), "Read complete independent summary");}
		for(Counts c:counts.values()){if(c!=null){require(c.checked[0] && c.checked[1] && c.n[0][0]==c.n[1][0] && c.n[1][1]-c.n[0][1]==c.gained-c.lost,"Both-arm reference conservation and net=gains-losses required");}}
	}
	public void old(String path){
		require(oldCases.isEmpty(),"One earlier55-call evidence table");final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		try{p.set(in.nextLine());require(p.terms()==32 && p.termEquals("case_id",0) && p.termEquals("origin",27),"Expected32-column mapped-call evidence");
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==32 && ids.add(p.parseString(0)),"Earlier case IDs must be unique");
				final Old c=new Old(p);final String key=c.assembly+'\t'+c.contig;ArrayList<Old> list=oldCases.get(key);
				if(list==null){list=new ArrayList<Old>();oldCases.put(key,list);}list.add(c);
			}require(ids.size()==55,"The declared historical characterization contains55 calls");
		}finally{require(!in.close(),"Read complete historical mapping evidence");}
	}
	public void delta(String path,String panel){
		require(panel.equals("P2") || panel.equals("P3"),"Explicit selected-panel GFF delta binding");
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{RrnaResourceIO.header(in,"change\tseqid\tsource\ttype\tstart\tstop\tscore\tstrand\tphase\tattributes",path);
			for(byte[] line;(line=in.nextLine())!=null;){p.set(line);require(p.terms()==10 && (p.termEquals("ADDED",0) || p.termEquals("REMOVED",0)),"Full original GFF9 delta required");
				final Feature f=new Feature(p);final String key=panel+'\t'+f.contig;ArrayList<Feature> list=changes.get(key);
				if(list==null){list=new ArrayList<Feature>();changes.put(key,list);}list.add(f);
			}
		}finally{require(!in.close(),"Read complete saved GFF changes");}
	}
	public void write(String prefix){
		require(!counts.isEmpty(),"Only reconciled evidence may be reported");RrnaResourceIO.fresh(prefix,".summary.tsv",".losses.tsv",".neighbors.tsv",".unmatched.tsv",".old_matches.tsv");
		final ByteStreamWriter sum=RrnaResourceIO.writer(prefix+".summary.tsv"),lost=RrnaResourceIO.writer(prefix+".losses.tsv"),near=RrnaResourceIO.writer(prefix+".neighbors.tsv"),fp=RrnaResourceIO.writer(prefix+".unmatched.tsv"),matches=RrnaResourceIO.writer(prefix+".old_matches.tsv");
		final ByteBuilder b=new ByteBuilder();
		try{
			lost.println("panel\tcompartment\tfamily\tassembly\tcontig\treference_id\ttruth_start1\ttruth_stop1\tstrand\tbaseline_call_id\tcall_start1\tcall_stop1\tcandidate_same_coordinate\tbaseline_GFF_attributes\tverifier_identity\tcausal_evidence");
			near.println("reference_id\tpanel\tcontig\tchange\ttype\tstart1\tstop1\tstrand\tGFF_score\ttruth_overlap_bp\tbaseline_call_overlap_bp\tbaseline_call_distance\tmodel\ttier\tattributes\tverifier_identity\tcausal_evidence");
			losses.sort((a,z)->a.refKey().compareTo(z.refKey()));
			for(Row r:losses){final Row call=calls.get(r.partner);require(call!=null && call.outcome.equals("TP") && call.partner.equals(r.id) && call.panel.equals(r.panel) && call.assembly.equals(r.assembly) && call.contig.equals(r.contig),"Lost reference must bind its exact baseline call partner");
				String attrs="NA";final ArrayList<Feature> list=changes.get(r.panel+'\t'+r.contig);
				if(list!=null){for(Feature f:list){if(f.change.equals("REMOVED") && f.same(call)){require(attrs.equals("NA"),"Ambiguous removed row for lost matched call");attrs=f.attrs;}
					final long distance=Math.max(0L,Math.max((long)call.start-f.stop,(long)f.start-call.stop));if(distance>1000){continue;}
					near.print(b.clear().append(r.id).tab().append(r.panel).tab().append(r.contig).tab().append(f.change).tab().append(f.type).tab().append(f.start).tab().append(f.stop).tab().append(f.strand).tab().append(f.score)
						.tab().append(overlap(r.start,r.stop,f.start,f.stop)).tab().append(overlap(call.start,call.stop,f.start,f.stop)).tab().append(distance).tab().append(f.model).tab().append(f.tier()).tab().append(f.attrs).append("\tNA_NOT_RECORDED\tSPATIAL_ONLY_DP_UNOBSERVED\n"));
				}}
				lost.print(b.clear().append(r.group()).tab().append(r.assembly).tab().append(r.contig).tab().append(r.id).tab().append(r.start).tab().append(r.stop).tab().append(r.strand).tab().append(call.id).tab().append(call.start).tab().append(call.stop)
					.tab().append(candidate.contains(call.coordinate())).tab().append(attrs).append("\tNA_NOT_RECORDED\tSPATIAL_ONLY_DP_UNOBSERVED\n"));
			}
			fp.println("panel\tcompartment\tfamily\tassembly\tcontig\tcall_id\tstart1\tstop1\tstrand\tbaseline_exact_call\tbaseline_exact_unmatched\told_relation\told_contig_cases\told_overlapping_cases\told_exact_cases\tinterpretation");
			matches.println("new_call_id\told_case_id\trelation\toverlap_bp\tsame_strand\told_origin\told_start1\told_stop1\told_primary_identity\tinterpretation");
			unmatched.sort((a,z)->a.coordinate().compareTo(z.coordinate()));
			for(Row r:unmatched){final Counts c=counts.get(r.group());if(oldFp.contains(r.coordinate())){c.retainedFp++;continue;}
				c.newFp++;final ArrayList<Old> old=oldCases.get(r.assembly+'\t'+r.contig);int same=0,overlap=0,exact=0;
				if(old!=null){for(Old o:old){same++;final int n=overlap(r.start,r.stop,o.start,o.stop);final boolean strand=r.strand.equals(o.strand),eq=strand && r.start==o.start && r.stop==o.stop;
					if(n>0 && strand){overlap++;}if(eq){exact++;}
					matches.print(b.clear().append(r.id).tab().append(o.id).tab().append(eq?"EXACT_COORDINATE":n>0 && strand?"SAME_STRAND_OVERLAP":"SAME_CONTIG_ONLY").tab().append(n).tab().append(strand).tab().append(o.origin).tab().append(o.start).tab().append(o.stop).tab().append(o.identity).append("\tPRIOR_CASE_EVIDENCE_NOT_NEW_CALL_IDENTITY_OR_FUNCTION\n"));
				}}
				final int relation=exact>0?0:overlap>0?1:same>0?2:3;c.relation[relation]++;
				fp.print(b.clear().append(r.group()).tab().append(r.assembly).tab().append(r.contig).tab().append(r.id).tab().append(r.start).tab().append(r.stop).tab().append(r.strand)
					.tab().append(baseline.contains(r.coordinate())).tab().append(false).tab().append(RELATIONS[relation]).tab().append(same).tab().append(overlap).tab().append(exact).append("\tREFERENCE_UNMATCHED_NOT_BIOLOGICAL_FP\n"));
			}
			sum.println("panel\tcompartment\tfamily\treferences\tbaseline_TP\tcandidate_TP\tlost_TP\tgained_TP\tbaseline_apparentFP\tcandidate_apparentFP\tnew_unmatched_coordinates\tretained_unmatched_coordinates\told_exact\told_overlap\told_contig_only\tnew_to_old55_contigs");
			for(String key:sorted(counts)){final Counts c=counts.get(key);require(c.checked[0] && c.checked[1],"Both saved summaries must validate every reported group");
				b.clear().append(key).tab().append(c.n[0][0]).tab().append(c.n[0][1]).tab().append(c.n[1][1]).tab().append(c.lost).tab().append(c.gained).tab().append(c.n[0][2]).tab().append(c.n[1][2]);
				if(key.endsWith("\tSSU")){require(c.newFp+c.retainedFp==c.n[1][2],"New+retained must conserve all candidate unmatched calls");b.tab().append(c.newFp).tab().append(c.retainedFp);for(long n:c.relation){b.tab().append(n);}}
				else{b.append("\tNA\tNA\tNA\tNA\tNA\tNA");}sum.print(b.nl());
			}
		}finally{RrnaResourceIO.finish(sum,lost,near,fp,matches);}
		System.err.println("MITO_SAVED_DELTA_PASS lost_references="+losses.size()+" no_alignment_no_rematching=true");
	}
	public static int overlap(int a,int z,int b,int y){require(a>0 && z>=a && b>0 && y>=b,"Inclusive genomic intervals required");return Math.max(0,Math.min(z,y)-Math.max(a,b)+1);}
	private static ArrayList<String> sorted(ObjectMap<String,?> m){final ArrayList<String> keys=new ArrayList<String>();for(String k:m.keys()){if(k!=null){keys.add(k);}}Collections.sort(keys);return keys;}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
	private static final class Row{
		Row(LineParser1 p){panel=p.parseString(1);assembly=p.parseString(2);contig=p.parseString(3);comp=p.parseString(4);family=p.parseString(5);id=p.parseString(7);start=p.parseInt(9);stop=p.parseInt(10);strand=p.parseString(11);outcome=p.parseString(12);partner=p.parseString(13);require(start>0 && stop>=start && (strand.equals("+") || strand.equals("-")),"Native audit coordinates must be valid");}
		String group(){return panel+'\t'+comp+'\t'+family;}
		String coordinate(){return panel+'\t'+assembly+'\t'+contig+'\t'+family+'\t'+start+'\t'+stop+'\t'+strand;}
		String refKey(){return group()+'\t'+assembly+'\t'+contig+'\t'+id;}
		boolean fp(){return outcome.equals("APPARENT_FP") || outcome.equals("DUPLICATE_APPARENT_FP");}
		final String panel,assembly,contig,comp,family,id,strand,outcome,partner;final int start,stop;
	}
	private static final class Old{
		Old(LineParser1 p){id=p.parseString(0);assembly=p.parseString(3);contig=p.parseString(4);start=p.parseInt(5);stop=p.parseInt(6);strand=p.parseString(7);origin=p.parseString(27);identity=p.parseString(18);require(start>0 && stop>=start,"Earlier characterized interval must be valid");}
		final String id,assembly,contig,strand,origin,identity;final int start,stop;
	}
	private static final class Feature{
		Feature(LineParser1 p){change=p.parseString(0);contig=p.parseString(1);type=p.parseString(3);start=p.parseInt(4);stop=p.parseInt(5);score=p.parseString(6);strand=p.parseString(7);attrs=p.parseString(9);final String m=PanelAnnotationGrade.colonField(attrs,"model");model=m==null?"NA":m;require(start>0 && stop>=start,"GFF delta retains valid original coordinates");}
		boolean same(Row r){return contig.equals(r.contig) && start==r.start && stop==r.stop && strand.equals(r.strand) && ((r.family.equals("SSU") && type.equals("rRNA")) || (r.family.equals("tRNA") && type.equals("tRNA")));}
		String tier(){return model.startsWith("mitochondria_coverage_")?"MITO_16S":type.equals("tRNA")?"TRNA":model.startsWith("euk18S")?"EUK18S":"OTHER";}
		final String change,contig,type,score,strand,attrs,model;final int start,stop;
	}
	private static final class Counts{final long[][] n=new long[2][3];final boolean[] checked=new boolean[2];long lost,gained,newFp,retainedFp;final long[] relation=new long[4];}
	private final ObjectMap<String,Row[]> refs=new ObjectMap<String,Row[]>(String.class,Row[].class);
	private final ObjectMap<String,Row> calls=new ObjectMap<String,Row>(String.class,Row.class);
	private final ObjectMap<String,Counts> counts=new ObjectMap<String,Counts>(String.class,Counts.class);
	private final ObjectMap<String,ArrayList<Old>> oldCases=new ObjectMap(String.class,ArrayList.class);
	private final ObjectMap<String,ArrayList<Feature>> changes=new ObjectMap(String.class,ArrayList.class);
	private final ObjectSet<String> baseline=new ObjectSet<String>(String.class),candidate=new ObjectSet<String>(String.class),oldFp=new ObjectSet<String>(String.class);
	private final ArrayList<Row> losses=new ArrayList<Row>(),unmatched=new ArrayList<Row>();
	private static final String[] REQUIRED={"audit","summary","old","p2delta","p3delta","outprefix"};
	private static final String KNOWN="|audit|summary|old|p2delta|p3delta|outprefix|";
	private static final String[] RELATIONS={"OLD_EXACT_COORDINATE","OLD_SAME_STRAND_OVERLAP","OLD_CONTIG_ONLY","NEW_TO_OLD55_CONTIG_SET"};
}
