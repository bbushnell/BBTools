package test.prok;

import java.io.File;
import java.util.Arrays;
import map.ObjectMap;
import prok.EukPanelGffValidation;
import prok.NcrnaMitoCallContexts;
import prok.RrnaResourceIO;

/** Merging, clipping, strict whole-contig boundary and source-orientation fixtures.
 * @author Brian, Ganyu */
public final class NcrnaMitoCallContextsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be fresh");
		final String small=NcrnaMitoCrossLibraryTest.sequence(199,7),large=NcrnaMitoCrossLibraryTest.sequence(2000,11),edge=NcrnaMitoCrossLibraryTest.sequence(200,13);
		write(root,"source.fa",">S assembly=A original_header=S\n"+small+"\n>L assembly=A original_header=L\n"+large+"\n>E assembly=A original_header=E\n"+edge+"\n");
		final String firstPins=NcrnaMitoCrossLibraryTest.bind(root,sha,new String[]{"source.fa"});final String hash=EukPanelGffValidation.hashes(firstPins).get(root+"/source.fa");
		final String mh="assembly\tcontig\tlength\tACGT\tGC\tGC_fraction_ACGT\texported_id\tsource_fasta\tsource_sha80\toriginal_header\toriginal_roster_sha80\tsource_hash_binding\n";
		write(root,"contigs.tsv",mh+source(root,hash,"S",199)+source(root,hash,"L",2000)+source(root,hash,"E",200));
		final String[] ch=new String[42];Arrays.fill(ch,"unused");ch[0]="case_id";ch[29]="audit_id";final String header=String.join("\t",ch)+"\n";
		final String data=call(root,hash,"C1","S",199,50,60,"+")+call(root,hash,"C2","L",2000,100,120,"-")+call(root,hash,"C3","L",2000,145,160,"+")
			+call(root,hash,"C4","L",2000,1000,1010,"+")+call(root,hash,"C5","E",200,195,200,"-");
		write(root,"cases.tsv",header+data);final String pins=NcrnaMitoCrossLibraryTest.bind(root,sha,new String[]{"source.fa","contigs.tsv","cases.tsv"});
		run(root,pins,"contexts","cases.tsv");
		final String summary=read(root+"/contexts.summary.tsv");check(summary.contains("5\t3\t2399\t4\t377\t1\t20\t200\tFORWARD_SOURCE\t0"),"Merged neighborhoods must conserve five calls with distinct original and measured denominators");
		final ObjectMap<String,String> seqs=new ObjectMap<String,String>(String.class,String.class);
		RrnaResourceIO.stream(root+"/contexts.fa.gz",r->{check(seqs.put(RrnaResourceIO.first(r.id),new String(r.bases,java.nio.charset.StandardCharsets.US_ASCII))==null,"Context IDs must be unique");});
		check(seqs.size()==4 && seqs.get("MCTX1").equals(edge.substring(174,200)) && seqs.get("MCTX2").equals(large.substring(79,180))
			&& seqs.get("MCTX3").equals(large.substring(979,1030)) && seqs.get("MCTX4").equals(small),"One-pass extraction must merge overlaps, clip ends, retain disjoint windows and never reverse-complement minus calls");
		final String memberships=read(root+"/contexts.calls.tsv");check(memberships.split("\n").length==6 && memberships.contains("C2\tMCTX2\tL\tA\tL\t100\t120\t-\t20\t40\t")
			&& memberships.contains("C3\tMCTX2\tL\tA\tL\t145\t160\t+\t65\t80\t"),"Merged calls must retain exact original and context-relative coordinates");
		check(read(root+"/contexts.contexts.tsv").contains("MCTX1\tE\tA\tE\t174\t199\t26\t26\t200\t200\tCALL_NEIGHBORHOOD"),"Whole-contig shortcut is strict less-than, never less-than-or-equal");
		write(root,"bad.tsv",header+data+call(root,hash,"C5","E",200,195,200,"-"));final String badPins=NcrnaMitoCrossLibraryTest.bind(root,sha,new String[]{"source.fa","contigs.tsv","bad.tsv"});
		boolean failed=false;try{run(root,badPins,"bad","bad.tsv");}catch(IllegalArgumentException e){failed=true;}check(failed && !new File(root+"/bad.summary.tsv").exists(),"Duplicate case IDs must fail without a completed summary");
		System.out.println("MITO_CALL_CONTEXTS_TEST_PASS merging disjoint_windows edge_clipping whole_boundary forward_orientation case_conservation duplicate_rejection");
	}
	private static void run(String root,String pins,String out,String cases){NcrnaMitoCallContexts.main(new String[]{"in="+root+"/source.fa","contigs="+root+"/contigs.tsv","cases="+root+"/"+cases,"verified="+pins,"outprefix="+root+"/"+out,"expected=5","expectedcontigs=3","flank=20","wholebelow=200"});}
	private static String source(String root,String hash,String id,int len){return "A\t"+id+'\t'+len+'\t'+len+"\t0\t0\t"+id+'\t'+root+"/source.fa\t"+hash+'\t'+id+'\t'+hash+"\tfixture\n";}
	private static String call(String root,String hash,String id,String source,int len,int start,int stop,String strand){final String[] p=new String[42];Arrays.fill(p,"NA");p[0]=id;p[1]="NUCLEAR_NEW_UNMATCHED";p[2]="A";p[3]=source;p[4]=""+start;p[5]=""+stop;p[6]=strand;p[7]=""+(stop-start+1);p[8]="mitochondria_insect";p[12]=""+len;p[13]=""+len;p[26]=root+"/source.fa";p[27]=hash;p[29]="fixture:"+id;return String.join("\t",p)+"\n";}
	private static void write(String root,String name,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+name,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
