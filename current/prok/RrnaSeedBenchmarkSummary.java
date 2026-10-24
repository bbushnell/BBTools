package prok;

import java.util.ArrayList;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import structures.ByteBuilder;

/** Reduce completed per-assembly grades without averaging precision or recall.
 * The euk5S combined row adds prok-control calls as FP, but no prok truth/FNs.
 * Seed coverage on prok controls stays diagnostic and is not pooled as euk truth.
 * @author Raiden
 */
public final class RrnaSeedBenchmarkSummary {
	public static void main(String[] args){
		String runs=null,out=null;int expected=-1;
		for(String arg:args){final int e=arg.indexOf('=');require(e>0,"Expected key=value");final String k=arg.substring(0,e),v=arg.substring(e+1);if(k.equals("runs")){runs=v;}else if(k.equals("out")){out=v;}else if(k.equals("expectedrows")){expected=Integer.parseInt(v);}else{throw new IllegalArgumentException(k);}}
		require(runs!=null&&out!=null&&expected>0,"Required runs=<comma-separated manifests> out= expectedrows=");new RrnaSeedBenchmarkSummary().run(runs,out,expected);
	}
	public void run(String paths,String out,int expected){
		assert(groups.size()==0):"One frozen collection per reducer";
		final ObjectSet<String> seen=new ObjectSet<String>(String.class);int rows=0;final LineParser1 p=new LineParser1('\t');
		for(String path:paths.split(",")){final ByteFile input=RrnaSourceTargetPrep.open(path);
			try{final byte[] h=input.nextLine();require(h!=null&&new String(h,java.nio.charset.StandardCharsets.UTF_8).equals("arm\tclade\tassembly\tgrade\tstages\tcallerSeconds"),"Unexpected completed-run manifest schema");
				for(byte[] line;(line=input.nextLine())!=null;){p.set(line);require(p.terms()==6,"Expected six completed-run columns");final String arm=p.parseString(0),clade=p.parseString(1),assembly=p.parseString(2);final String grade=p.parseString(3),stages=p.parseString(4);final double seconds=p.parseDouble(5);
					require(Double.isFinite(seconds)&&seconds>=0,"Caller time must be a finite nonnegative observation");final Row row=loadGrade(grade,arm,clade,assembly);require(seen.add(arm+"\t"+clade+"\t"+assembly+"\t"+row.family),"Duplicate completed genome/arm/family");
					if(!stages.equals("NA")){loadStages(stages,row,arm,clade,assembly);}else{require(row.family.equals("prok5S"),"Every euk5S arm requires measured stages");}
					add(arm,clade,row,seconds,false);add(arm,"combined",row,seconds,true);rows++;
				}
			}finally{require(!input.close(),"Completed-run manifest read failed");}
		}
		require(rows==expected,"Completed-run row count "+rows+" differs from declared "+expected);
		final ByteStreamWriter output=RrnaSourceTargetPrep.writer(out);final ByteBuilder b=new ByteBuilder();
		try{output.println("arm\tclade\tfamily\tgenomes\tTP\tFP\tFN\tprecision\trecall\texact5\texact3\tbothExact\texact5Fraction\texact3Fraction\tbothExactFraction\tcontrolOverlapAnyStrand\tgenomeBases\tstageReferenceLoci\tseedCoveredReferenceLoci\tseedCoverage\twindowEvents\tuniqueWindows\twindowsPerMb\talignments\tcallerSeconds\tduplicateCalls");
			for(Group g:ordered){b.clear().append(g.arm).tab().append(g.clade).tab().append(g.family).tab().append(g.genomes).tab().append(g.tp).tab().append(g.fp).tab().append(g.fn).tab();RrnaSeedBenchmarkGrade.ratio(b,g.tp,g.tp+g.fp);b.tab();RrnaSeedBenchmarkGrade.ratio(b,g.tp,g.tp+g.fn);b.tab().append(g.exact5).tab().append(g.exact3).tab().append(g.both).tab();RrnaSeedBenchmarkGrade.ratio(b,g.exact5,g.tp);b.tab();RrnaSeedBenchmarkGrade.ratio(b,g.exact3,g.tp);b.tab();RrnaSeedBenchmarkGrade.ratio(b,g.both,g.tp);b.tab().append(g.overlap).tab();
				if(g.stageGenomes==0){b.append("NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA");}else{
					require(g.stageGenomes==g.genomes,"No partial-stage aggregate is allowed");b.append(g.bases).tab().append(g.reference).tab().append(g.covered).tab();RrnaSeedBenchmarkGrade.ratio(b,g.covered,g.reference);b.tab().append(g.events).tab().append(g.windows).tab();if(g.bases>0){b.append(g.windows*1_000_000.0/g.bases,8);}else{b.append("NA");}b.tab().append(g.alignments);
				}b.tab().append(g.seconds,6).tab().append(g.duplicates).nl();output.print(b);
			}
		}finally{RrnaCmWindowParity.finish(output);}
		System.err.println("SEED_BENCHMARK_SUMMARY_PASS rows="+rows+" groups="+ordered.size()+"; pooled ratios from counts, prok euk5S calls are FP only");
	}
	public static Row loadGrade(String path,String arm,String clade,String assembly){
		final LineParser1 p=oneRow(path,22);identity(p,arm,clade,assembly);final Row r=new Row();r.family=p.parseString(3);r.negative=p.parseInt(4)==1;
		require(r.family.equals("euk5S") || r.family.equals("prok5S") || r.family.equals("legacy5S"),"Unknown graded family");r.tp=p.parseLong(8);r.fp=p.parseLong(9);r.fn=p.parseLong(10);r.exact5=p.parseLong(13);r.exact3=p.parseLong(14);r.both=p.parseLong(15);r.duplicates=p.parseLong(19);r.overlap=p.parseLong(20);
		require(r.tp>=0&&r.fp>=0&&r.fn>=0&&r.tp+r.fp==p.parseLong(7)&&r.tp+r.fn==p.parseLong(5),"Grade conservation failure");
		require(r.exact5>=r.both&&r.exact3>=r.both&&r.both>=0&&r.exact5<=r.tp&&r.exact3<=r.tp,"Endpoint counts must be TP subsets");
		require(!r.negative||(r.family.equals("euk5S")&&r.tp==0&&r.fn==0),"Negative controls cannot contribute euk truth");return r;
	}
	static void loadStages(String path,Row r,String arm,String clade,String assembly){
		assert(r!=null):"Stage metrics are joined to an existing final-call grade";final LineParser1 p=oneRow(path,16);identity(p,arm,clade,assembly);r.stages=true;r.bases=p.parseLong(3);r.reference=p.parseLong(5);r.covered=p.parseLong(6);r.events=p.parseLong(9);r.windows=p.parseLong(10);r.alignments=p.parseLong(15);
		require(r.bases>=0&&r.covered>=0&&r.covered<=r.reference&&r.windows>=0&&r.windows<=r.events&&r.alignments>=0,"Invalid measured-stage counts");
		require(r.negative||r.reference==r.tp+r.fn,"Seed coverage and final-call grade must use the same biological truth");
	}
	static LineParser1 oneRow(String path,int columns){
		final ByteFile in=RrnaSourceTargetPrep.open(path);final LineParser1 p=new LineParser1('\t');
		try{final byte[] header=in.nextLine(),row=in.nextLine();require(header!=null&&row!=null&&in.nextLine()==null,"Expected one summary header and one data row: "+path);p.set(header);require(p.terms()==columns&&p.termEquals("arm",0)&&p.termEquals("clade",1)&&p.termEquals("assembly",2),"Summary schema changed: "+path);p.set(row);require(p.terms()==columns,"Summary column count changed");return p;
		}finally{require(!in.close(),"Summary read failed");}
	}
	static void identity(LineParser1 p,String arm,String clade,String assembly){require(p.termEquals(arm,0)&&p.termEquals(clade,1)&&p.termEquals(assembly,2),"Summary belongs to a different arm/clade/assembly");}
	void add(String arm,String clade,Row r,double seconds,boolean pooled){
		assert(r!=null):"Only completed per-assembly rows may be aggregated";final String key=arm+"\t"+clade+"\t"+r.family;Group g=groups.get(key);if(g==null){g=new Group(arm,clade,r.family);groups.put(key,g);ordered.add(g);}
		g.genomes++;g.tp+=r.tp;g.fp+=r.fp;g.fn+=r.fn;g.exact5+=r.exact5;g.exact3+=r.exact3;g.both+=r.both;g.overlap+=r.overlap;g.duplicates+=r.duplicates;g.seconds+=seconds;
		if(r.stages){g.stageGenomes++;g.bases+=r.bases;g.events+=r.events;g.windows+=r.windows;g.alignments+=r.alignments;if(!pooled||!r.negative){g.reference+=r.reference;g.covered+=r.covered;}}
	}
	static void require(boolean b,String why){RrnaSeedBenchmarkGrade.require(b,why);}
	public static class Row{String family;public boolean negative,stages;public long tp,fp,fn,exact5,exact3,both,overlap,duplicates,bases,reference,covered,events,windows,alignments;}
	static final class Group extends Row{Group(String a,String c,String f){arm=a;clade=c;family=f;}final String arm,clade;long genomes,stageGenomes;double seconds;}
	final ObjectMap<String,Group> groups=new ObjectMap<String,Group>(String.class,Group.class);final ArrayList<Group> ordered=new ArrayList<Group>();
}
