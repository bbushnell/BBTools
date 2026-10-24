package test.prok;

import prok.*;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;

/** Geometry, denominator and full-member preservation checks for Q3.
 * @author Raiden */
public final class EukPanelLociTest {
	public static void main(String[] args){
		require(args.length==1 && args[0].startsWith("outdir="), "Expected fresh outdir");
		final String root=args[0].substring(7);require(new File(root).mkdir(), "Fixture directory must be fresh");
		geometry();
		final String a=row(1, "RF00177", "SSU", 100, 199, 200, true), b=partial(row(2, "RF01960", "SSU", 151, 199, 210, true));
		final String c=row(3, "RF00177", "SSU", 300, 399, 20, true), d=row(4, "RF00177", "SSU", 990, 1010, 60, false);
		write(root+"/hits.tsv", EukPanelLoci.HEADER+'\n'+a+'\n'+b+'\n'+c+'\n'+d+'\n');
		write(root+"/models.tsv", "bundle_sha80\taccession\tfamily\tsubtype\tclen\tga\tsource_cm\tsource_cm_sha80\n"+
			model("RF00177","SSU")+model("RF01960","SSU")+model("RF00005","tRNA"));
		write(root+"/scans.tsv", "selection_bin\t"+EukPanelCmProposal.SUMMARY_HEADER+'\n'+scan("RF00177",3)+scan("RF01960",1)+scan("RF00005",0));
		run(root, "hits.tsv", "one", 4);verify(root, "one", new String[]{a,b,c}, d);
		write(root+"/reversed.tsv", EukPanelLoci.HEADER+'\n'+d+'\n'+c+'\n'+b+'\n'+a+'\n');
		run(root, "reversed.tsv", "two", 4);
		for(String suffix:new String[]{".loci.tsv.gz",".locus_counts.tsv",".member_distribution.tsv",".locus_cases.tsv"}){
			require(read(root+"/one"+suffix).equals(read(root+"/two"+suffix)), "Reordering hits must preserve locus IDs and summaries: "+suffix);
		}
		write(root+"/duplicate.tsv", EukPanelLoci.HEADER+'\n'+a+'\n'+a+'\n');expect(root,"duplicate.tsv","bad_duplicate",2,"Duplicate run/line");
		write(root+"/missing.tsv", EukPanelLoci.HEADER+'\n'+a+'\n');expect(root,"missing.tsv","bad_count",1,"scan denominator");
		System.out.println("PASS tandem bridge deferral, residual spanning flag, member-half threshold, full tie-break, input permutation, preserved raw members/invalid rows, zero-hit counts, GA strata and fail-loud duplicate/denominator checks");
	}
	static void geometry(){
		final EukPanelLoci.Hit left=hit(1,100,199,100), right=hit(2,300,399,90), bridge=hit(3,100,399,500);
		ArrayList<EukPanelLoci.Locus> loci=build(left,right,bridge);
		require(loci.size()==3 && bridge.spanning && left.locus!=right.locus && bridge.locus.seed==bridge, "High-score spanning seed must not absorb two separate copies");
		final EukPanelLoci.Hit a=hit(1,100,199,100), b=hit(2,180,279,90), wide=hit(3,100,279,500);
		loci=build(wide,b,a);require(loci.size()==2 && wide.spanning && wide.locus==a.locus, "Deferred spanning hit joins earliest qualifying non-spanning seed");
		final EukPanelLoci.Hit p=hit(1,1,101,50), q=hit(2,52,152,40), r=hit(3,51,151,30);
		require(!EukPanelLoci.qualifies(q,p) && EukPanelLoci.qualifies(r,p), "Odd length101 requires51bp, not50bp");
		final EukPanelLoci.Hit x=hit(10,1,100,50), y=hit(9,1,100,50);
		loci=build(x,y);require(loci.size()==1 && loci.get(0).seed==y, "Final source-line tie-break must be deterministic");
		final EukPanelLoci.Hit s=hit(1,100,199,90), tiny=hit(2,120,149,20);
		loci=build(tiny,s);require(loci.size()==1 && tiny.locus==s.locus, "Membership is half the member, not reciprocal half");
		final EukPanelLoci.Hit seed=hit(1,100,299,90), u=hit(2,100,149,80), v=hit(3,250,299,80);
		final EukPanelLoci.Locus forced=new EukPanelLoci.Locus(seed,1);forced.add(u);forced.add(v);EukPanelLoci.findAmbiguity(forced);
		require(forced.witnessA!=null, "Two separated absorbed members need explicit witnesses");
		final EukPanelLoci.Hit container=hit(1,100,299,500), belowA=hit(2,100,149,20), belowB=hit(3,250,299,20);
		loci=build(belowA,container,belowB);require(loci.size()==1 && !container.spanning, "Below-GA pieces do not establish independent tandem copies");
		final EukPanelLoci.Hit full=hit(1,100,299,500), partialA=partialHit(2,100,149), partialB=partialHit(3,250,299);
		loci=build(partialA,full,partialB);require(loci.size()==1 && !full.spanning, "GA-passing partial model pieces must not demote a containing full hit");
		final EukPanelLoci.Hit lowerFull=hit(1,100,299,60), higherPartial=partialHit(2,250,299);
		loci=build(higherPartial,lowerFull);require(loci.size()==1 && loci.get(0).seed==lowerFull && EukPanelLoci.subtype(loci.get(0).members)==higherPartial,"Full GA geometry seed absorbs a stronger partial without overriding subtype evidence");
		final EukPanelLoci.Hit below=hit(3,250,299,40);
		require(EukPanelLoci.RANK.compare(lowerFull,higherPartial)<0 && EukPanelLoci.RANK.compare(higherPartial,below)<0,"Seed evidence classes order full GA, partial GA, below GA");
		loci=build(below);require(EukPanelLoci.subtype(loci.get(0).members)==null,"Below-GA-only loci must not acquire supported subtype labels");
	}
	static String partial(String row){return row.replace("\t1\t100\t100\t100\t1.0\t","\t10\t90\t100\t81\t0.81\t").replace("\ttrue\t0\t0\t","\tfalse\t9\t10\t");}
	static EukPanelLoci.Hit partialHit(int line,long start,long end){
		final String row=partial(row(line,"RF00177","SSU",start,end,100,true));
		final byte[] raw=row.getBytes(StandardCharsets.UTF_8);final LineParser1 p=new LineParser1('\t');p.set(raw);
		final EukPanelLoci.Hit h=new EukPanelLoci.Hit(raw,p);require(h.h.gaPass && !h.h.full,"Substantial-witness negative fixture must really be partial");return h;
	}
	static ArrayList<EukPanelLoci.Locus> build(EukPanelLoci.Hit... hits){
		final ArrayList<EukPanelLoci.Hit> group=new ArrayList<EukPanelLoci.Hit>();Collections.addAll(group,hits);
		final ArrayList<EukPanelLoci.Locus> out=new ArrayList<EukPanelLoci.Locus>();EukPanelLoci.buildPartition(group,out);return out;
	}
	static EukPanelLoci.Hit hit(int line,long start,long end,double bits){final LineParser1 p=new LineParser1('\t');final byte[] raw=row(line,"RF00177","SSU",start,end,bits,true).getBytes(StandardCharsets.UTF_8);p.set(raw);return new EukPanelLoci.Hit(raw,p);}
	static String row(int line,String accession,String family,long start,long end,double bits,boolean valid){
		return "A\trun\t/raw.tbl\t00000000000000000000\t"+line+"\t00000000000000000000\t"+accession+'\t'+family+"\ttest\tchr\tchr\tUNKNOWN\tEXACT\t"+(valid?"VALID_INTERVAL":"OUT_OF_RANGE")+
			"\t1\t100\t100\t100\t1.0\t"+start+'\t'+end+"\t+\tno\t!\t"+bits+"\t1e-10\t50\t"+(bits-50)+"\t0\ttest\t1000\t"+(valid?(start==1?"LEFT":"INTERNAL"):"UNASSESSABLE")+"\ttrue\t0\t0\t"+(bits>=50)+"\t0\t0\t0\tREVIEW_ONLY";
	}
	static String model(String a,String f){return "00000000000000000000\t"+a+'\t'+f+"\ttest\t100\t50\t/model.cm\t00000000000000000000\n";}
	static String scan(String a,int n){return "test\tA\trun\t"+a+"\t00000000000000000000\t"+n+"\t0\t0\t0\t0\t0\t0\t0\t0\n";}
	static void run(String root,String file,String out,int n){EukPanelLoci.main(new String[]{"in="+root+'/'+file,"models="+root+"/models.tsv","scans="+root+"/scans.tsv","outprefix="+root+'/'+out,"expectedhits="+n,"expectedmodels=3","expectedscans=3"});}
	static void verify(String root,String out,String[] raw,String invalid){
		final String prefix=root+'/'+out;final String[] members=read(prefix+".locus_members.tsv.gz").split("\n");
		require(members.length==4, "All three valid hits remain members");for(int i=0;i<raw.length;i++){require(members[i+1].startsWith(raw[i]+'\t'), "Every original field must survive");}
		require(read(prefix+".invalid_hits.tsv.gz").equals(EukPanelLoci.HEADER+'\n'+invalid+'\n'), "Invalid raw row retained separately");
		final String counts=read(prefix+".locus_counts.tsv");require(counts.contains("test\tSSU\t2\t3\t1\t1\t0\t0\n") && counts.contains("test\ttRNA\t0\t0\t0\t0\t0\t0\n"), "Counts include GA strata and zero-hit family");
		final String[] rows=read(prefix+".loci.tsv.gz").split("\n");
		require(rows[0].equals(EukPanelLoci.LOCUS_HEADER),"Geometry seed and supported subtype have explicit independent columns");
		final String[] first=rows[1].split("\t"), second=rows[2].split("\t");
		require(first.length==26 && first[9].equals("RF00177") && first[25].equals("RF01960") && first[4].equals("100") && first[5].equals("199"),"Written geometry uses full seed while subtype uses stronger partial member");
		require(second[25].equals("NA") && first[23].equals("DESCRIPTIVE_ONLY"),"No below-GA subtype or automatic truth admission");
	}
	static void expect(String root,String file,String out,int n,String message){boolean failed=false;try{run(root,file,out,n);}catch(IllegalArgumentException e){require(e.getMessage().contains(message), "Unexpected failure: "+e);failed=true;}require(failed,"Missing failure: "+message);}
	static void write(String path,String text){final ByteStreamWriter out=RrnaResourceIO.writer(path);out.print(text);RrnaResourceIO.finish(out);}
	static String read(String path){final ByteFile in=RrnaResourceIO.open(path);final StringBuilder b=new StringBuilder();try{for(byte[] line;(line=in.nextLine())!=null;){b.append(new String(line,StandardCharsets.UTF_8)).append('\n');}}finally{require(!in.close(),"Fixture read failed");}return b.toString();}
	static void require(boolean ok,String why){EukPanelLoci.require(ok,why);}
}
