package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import prok.EukPanelTruthFreeze;
import prok.NcrnaMitoCoverageSelect;

/** Exact-cache, fragment, alias, deterministic cover and matrix-binding fixtures.
 * @author Brian, Ganyu */
public final class NcrnaMitoCoverageSelectTest {
	public static void main(String[] args)throws Exception{
		check(args.length==2 && args[0].startsWith("outdir=") && args[1].startsWith("sha80="),"Require fresh outdir and canonical sha80 executable");
		final String root=args[0].substring(7),sha=args[1].substring(6);check(new File(root).mkdir() && new File(sha).canExecute(),"Fixture directory must be fresh");
		final String a=NcrnaMitoCrossLibraryTest.sequence(400,17),b=NcrnaMitoCrossLibraryTest.sequence(400,29);
		write(root,"corpus.fa",record("A1","uA",4,false,a)+record("A2","uA",999,true,a)+record("B1","uB",4,false,b));
		write(root,"ssu.fa",">A n=2\n"+a+"\n>B n=1\n"+b+"\n");write(root,"nohit.fa",">DUP n=1\n"+a+"\n");
		final String mh="arm\tlocus_id\tunique_id\tsource_taxid\tsource_length\tpartial\tpseudo\tsource_product\tlength_bin\tlineage_rank\tlineage\tuniversal_preference\tbest_ssu_cm\tcm_concordance\tassignment_status\tcluster_key\tassigned_model\tquantum_identity\ttrace_identity\tpaired_columns\tmember_paired_fraction\tmodel_paired_fraction\n";
		write(root,"ssu.members.tsv",mh+member("ssu","A1","uA",4,false,"A")+member("ssu","A2","uA",999,true,"A")+member("ssu","B1","uB",4,false,"B"));
		write(root,"nohit.members.tsv",mh+member("nohit","OTHER","uA",4,false,"DUP"));
		final String[] header=new String[42];Arrays.fill(header,"unused");header[1]="kind";header[29]="audit_id";header[37]="host_taxid";
		write(root,"cases.tsv",String.join("\t",header)+"\n"+panel("P1","C1",1)+"\n"+panel("P2","C2",2)+"\n");
		write(root,"panel.fa",">P1 assembly=P contig=C1 start=1 stop=80 strand=+\n"+a.substring(50,130)+"\n>P2 assembly=P contig=C2 start=1 stop=80 strand=+\n"+repeat('N',80)+"\n");
		write(root,"truth.tsv",EukPanelTruthFreeze.HEADER+"\n"+truth("1","C1")+"\n"+truth("2","C2")+"\n");NcrnaMitoCrossLibraryTest.tree(root+"/tree.tsv");
		final String pins=NcrnaMitoCrossLibraryTest.bind(root,sha,FILES);
		run(root,pins,"serial","full",null,1);run(root,pins,"parallel","full",null,2);
		for(String suffix:new String[]{".matrix.tsv.gz",".members.tsv.gz",".models.tsv",".selection.tsv",".selected.fa",".composition.tsv",".summary.tsv",".strata.tsv"}){
			check(read(root+"/serial"+suffix).equals(read(root+"/parallel"+suffix)),"Serial/parallel output must be identical: "+suffix);}
		final String[] result=read(root+"/serial.summary.tsv").split("\n")[1].split("\t");
		check(result[0].equals("5") && result[1].equals("4") && result[2].equals("3") && result[3].equals("2") && result[5].equals("3") && result[6].equals("1") && result[7].equals("4"),"Duplicate copies and duplicate candidate models must not inflate unique coverage; uncovered Ns remain");
		check(result[12].equals("full") && result[13].equals("2") && result[14].equals("8") && result[15].equals("12"),"Full-mode work counts must exclude cached model aliases while matrix rows retain them");
		check(read(root+"/serial.selected.fa").contains(">mitochondria_coverage_ssu_A") && !read(root+"/serial.selected.fa").contains("nohit_DUP"),"Frozen first-model tie keeps the canonical alias without discarding its provenance");
		check(read(root+"/serial.models.tsv").contains("nohit:DUP\tmitochondria_coverage_nohit_DUP\t400\tssu:A"),"Duplicate candidate aliases must remain explicit");
		check(read(root+"/serial.composition.tsv").contains("nohit:DUP\tSOURCE_CLUSTER_FOR_SELECTED_ssu:A"),"Aliased original source-cluster composition must not disappear");
		check(read(root+"/serial.members.tsv.gz").contains("PANEL_INCLUDED\tP2\t") && read(root+"/serial.members.tsv.gz").contains("UNCOVERED_AT_CUT"),"Every panel locus remains even when no candidate can cover it");
		check(read(root+"/serial.strata.tsv").contains("PANEL_INCLUDED:PARTIAL\t2\t2\t1\t1\t1\t1"),"Shorter-query fragment coverage must be reported separately, including residuals");
		final String matrix=root+"/serial.matrix.tsv.gz";final String[] extended=Arrays.copyOf(FILES,FILES.length+1);extended[FILES.length]="serial.matrix.tsv.gz";
		final String matrixPins=NcrnaMitoCrossLibraryTest.bind(root,sha,extended);run(root,matrixPins,"reselect","select",matrix,1);
		check(read(root+"/serial.selected.fa").equals(read(root+"/reselect.selected.fa")) && read(root+"/serial.members.tsv.gz").equals(read(root+"/reselect.members.tsv.gz")),"A bound reusable matrix must reproduce selection without realignment");
		final String[] reused=read(root+"/reselect.summary.tsv").split("\n")[1].split("\t");
		check(reused[12].equals("select") && reused[13].equals("2") && reused[14].equals("0") && reused[15].equals("12"),"Matrix-only selection must report zero new alignment calls, with all cached rows accounted for");
		run(root,pins,"pilot","pilot",null,1);check(!new File(root+"/pilot.selected.fa").exists(),"Timing pilot must not silently produce a resource selection");
		final String pilot=read(root+"/pilot.pilot.tsv");check(pilot.startsWith("unique_id\tlength\tconstruction_member\tpanel_member\tplanned_alignments\n") && pilot.split("\n").length==4,"Frozen pilot membership must label planned work and include both panel sequences plus one length-selected construction sequence");
		final String[] pilotFiles=Arrays.copyOf(FILES,FILES.length+1);pilotFiles[FILES.length]="pilot.matrix.tsv.gz";final String pilotPins=NcrnaMitoCrossLibraryTest.bind(root,sha,pilotFiles);
		boolean failed=false;try{run(root,pilotPins,"bad_pilot_select","select",root+"/pilot.matrix.tsv.gz",1);}catch(IllegalArgumentException e){failed=true;}check(failed,"A partial pilot matrix cannot substitute for full coverage evidence");
		check(Arrays.equals(NcrnaMitoCoverageSelect.greedy(new long[][]{{3},{6},{8}},5),new int[]{0,1,2}),"Greedy uses uncovered unique bits and frozen order for ties, leaving unreachable members");
		final long[][] many=new long[95][2];for(int i=0;i<95;i++){many[i][i>>>6]=1L<<(i&63);}check(NcrnaMitoCoverageSelect.greedy(many,95).length==95,"Expected10–20models must never become a forced cap; tailword bits conserve coverage");
		System.out.println("MITO_COVERAGE_SELECT_TEST_PASS copies unique_cache aliases fragments residuals deterministic matrix_reuse pilot_guard uncapped_bitset");
	}
	private static void run(String root,String pins,String name,String mode,String matrix,int threads){
		final ArrayList<String> args=new ArrayList<String>();for(int i=0;i<KEYS.length;i++){args.add(KEYS[i]+"="+root+"/"+FILES[i]);}
		args.add("verified="+pins);args.add("outprefix="+root+"/"+name);args.add("mode="+mode);args.add("expected=3");args.add("expectedpanel=2");args.add("expectedssu=2");args.add("expectednohit=1");args.add("expectedssucopies=3");args.add("expectednohitcopies=1");args.add("pilotsize=1");args.add("t="+threads);if(matrix!=null){args.add("matrix="+matrix);}NcrnaMitoCoverageSelect.main(args.toArray(new String[0]));
	}
	private static String member(String role,String id,String uid,int tid,boolean partial,String model){final String[] p=new String[22];Arrays.fill(p,"NA");p[0]=role;p[1]=id;p[2]=uid;p[3]=""+tid;p[4]="400";p[5]=""+partial;p[6]="false";p[14]="kept";p[16]=model+" n=1";return String.join("\t",p)+"\n";}
	private static String panel(String id,String contig,int locus){final String[] p=new String[42];Arrays.fill(p,"NA");p[0]=id;p[1]="ORGANELLE_SSU_TRUTH";p[2]="P";p[3]=contig;p[4]="1";p[5]="80";p[6]="+";p[29]="LOCUS:"+locus;p[37]="4";return String.join("\t",p);}
	private static String truth(String id,String contig){final String[] p=new String[35];Arrays.fill(p,"NA");p[0]=id;p[1]="P";p[2]=contig;p[3]="+";p[4]="1";p[5]="80";p[6]="SSU";p[15]="0";p[16]="0";p[27]="PARTIAL_INTERNAL";return String.join("\t",p);}
	private static String record(String id,String uid,int tid,boolean partial,String seq){return ">"+id+" unique="+uid+" tid|"+tid+" partial="+partial+" pseudo=false\n"+seq+"\n";}
	private static String repeat(char c,int n){final char[] a=new char[n];Arrays.fill(a,c);return new String(a);}
	private static void write(String root,String name,String text){NcrnaMitoCrossLibraryTest.write(root+"/"+name,text);}
	private static String read(String path){return NcrnaMitoCrossLibraryTest.read(path);}
	private static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
	private static final String[] KEYS={"corpus","ssumodels","nohitmodels","ssumembers","nohitmembers","panelcases","panelgenes","truth","tree"};
	private static final String[] FILES={"corpus.fa","ssu.fa","nohit.fa","ssu.members.tsv","nohit.members.tsv","cases.tsv","panel.fa","truth.tsv","tree.tsv"};
}
